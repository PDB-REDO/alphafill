/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2021 Maarten L. Hekkelman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// This code is originally written for mini-ibs, a content management system

#include "https-client.hpp"
#include "revision.hpp"

#include <boost/asio.hpp>
#include <boost/asio/ssl.hpp>
#include <boost/version.hpp>
#include <iostream>
#include <mcfp/mcfp.hpp>
#include <zeep/http/message-parser.hpp>
#include <zeep/streambuf.hpp>
#include <zeep/unicode-support.hpp>

// --------------------------------------------------------------------
// Sometimes we need to fetch media that is not available yet

using boost::asio::ip::tcp;

template <typename SocketType>
class client_base
{
  public:
	using socket_type = SocketType;

	virtual ~client_base() = default;

	[[nodiscard]] bool done() const { return m_done; }
	zeep::http::reply get_reply() { return m_reply_parser.get_reply(); }

  protected:
	virtual socket_type &get_socket() = 0;

	explicit client_base(const zeep::uri &url)
		: m_req({ "GET", url })
		, m_verbose(mcfp::config::instance().has("m_verbose"))
	{
	}

	void send_request()
	{
		std::vector<boost::asio::const_buffer> buffers;
		for (auto &buffer : m_req.to_buffers())
			buffers.emplace_back(buffer.data(), buffer.size());

		boost::asio::async_write(get_socket(),
			buffers,
			[this](const boost::system::error_code &error, std::size_t /*length*/)
			{
				if (not error)
					receive_response();
				else if (m_verbose)
					std::cout << "Write failed: " << error.message() << "\n";
			});
	}

	void receive_response()
	{
		boost::asio::async_read(get_socket(),
			boost::asio::buffer(m_buffer),
			[this](const boost::system::error_code &error, std::size_t length)
			{
				if (error and error != boost::asio::error::eof)
				{
					if (m_verbose)
						std::cout << "Read failed: " << error.message() << "\n";
					return;
				}

				zeep::char_streambuf sb(m_buffer.data(), length);

				auto r = m_reply_parser.parse(sb);
				if (r == true or error == boost::asio::error::eof)
					m_done = true;
				else
					receive_response();
			});
	}

	std::array<char, 4096> m_buffer{};
	const zeep::http::request m_req;
	bool m_done = false, m_verbose = false;
	zeep::http::reply_parser m_reply_parser;
};

class client : public client_base<tcp::socket>
{
  public:
	client(boost::asio::io_context &io_context,
		const tcp::resolver::results_type &endpoints,
		const std::string &req)
		: client_base(req)
		, m_socket(io_context)
	{
		connect(endpoints);
	}

  protected:
	socket_type &get_socket() override { return m_socket; }

  private:
	void connect(const tcp::resolver::results_type &endpoints)
	{
		boost::asio::async_connect(m_socket, endpoints,
			[this](const boost::system::error_code &error,
				const tcp::endpoint & /*endpoint*/)
			{
				if (not error)
					send_request();
				else if (m_verbose)
					std::cout << "Connect failed: " << error.message() << "\n";
			});
	}

	tcp::socket m_socket;
};

class ssl_client : public client_base<boost::asio::ssl::stream<tcp::socket>>
{
  public:
	ssl_client(boost::asio::io_context &io_context,
		boost::asio::ssl::context &context,
		const tcp::resolver::results_type &endpoints,
		const std::string &req)
		: client_base(req)
		, m_socket(io_context, context)
	{
		m_socket.set_verify_mode(boost::asio::ssl::verify_peer);
		m_socket.set_verify_callback(
			[this](auto &&preverified, auto &&ctx)
			{ return verify_certificate(
				  std::forward<decltype(preverified)>(preverified),
				  std::forward<decltype(ctx)>(ctx)); });

		connect(endpoints);
	}

  protected:
	socket_type &get_socket() override { return m_socket; }

  private:
	bool verify_certificate(bool preverified,
		boost::asio::ssl::verify_context & /*ctx*/)
	{
		// // The verify callback can be used to check whether the certificate that is
		// // being presented is valid for the peer. For example, RFC 2818 describes
		// // the steps involved in doing this for HTTPS. Consult the OpenSSL
		// // documentation for more details. Note that the callback is called once
		// // for each certificate in the certificate chain, starting from the root
		// // certificate authority.

		// // In this example we will simply print the certificate's subject name.
		// char subject_name[256];
		// X509 *cert = X509_STORE_CTX_get_current_cert(ctx.native_handle());
		// X509_NAME_oneline(X509_get_subject_name(cert), subject_name, 256);
		// std::cout << "Verifying " << subject_name << "\n";

		return preverified;
	}

	void connect(const tcp::resolver::results_type &endpoints)
	{
		boost::asio::async_connect(m_socket.lowest_layer(), endpoints,
			[this](const boost::system::error_code &error,
				const tcp::endpoint & /*endpoint*/)
			{
				if (not error)
					handshake();
				else if (m_verbose)
					std::cout << "Connect failed: " << error.message() << "\n";
			});
	}

	void handshake()
	{
		m_socket.async_handshake(boost::asio::ssl::stream_base::client,
			[this](const boost::system::error_code &error)
			{
				if (not error)
					send_request();
				else if (m_verbose)
					std::cout << "Handshake failed: " << error.message() << "\n";
			});
	}

	boost::asio::ssl::stream<tcp::socket> m_socket;
};

zeep::http::reply send_request(zeep::http::request &req, const zeep::uri &url)
{
	namespace ssl = boost::asio::ssl;
	using ssl_socket = ssl::stream<tcp::socket>;

	auto host = url.get_host();
	auto port = url.get_port();
	if (port == 0)
		port = url.get_scheme() == "http" ? 80 : 443;

	boost::asio::io_context io_context;
	tcp::resolver resolver(io_context);
	tcp::resolver::results_type endpoints = resolver.resolve(host, std::to_string(port));

	// prepare a request
	req.set_header("Host", std::format("{}:{}", host, port));

	if (req.get_header("accept").empty())
		req.set_header("Accept", "*/*");

	if (req.get_header("user-agent").empty())
		req.set_header("User-Agent", std::format("{}/{}", kProjectName, kVersionNumber));

	std::vector<boost::asio::const_buffer> req_buffer;
	for (auto &buffer : req.to_buffers())
		req_buffer.emplace_back(buffer.data(), buffer.size());

	auto reader = [&, is_head = zeep::iequals(req.get_method(), "HEAD")](auto &socket)
	{
		zeep::http::reply result;
		zeep::http::reply_parser p;

		for (;;)
		{
			std::array<char, 4096> buf{};
			boost::system::error_code error{};

			size_t len = socket.read_some(boost::asio::buffer(buf), error);

			zeep::char_streambuf sb(buf.data(), len);

			auto r = p.parse(sb);

			if (r == true or error == boost::asio::error::eof or len == 0 or (sb.in_avail() == 0 and is_head))
			{
				result = p.get_reply();
				break;
			}
			else if (error)
			{
				if (mcfp::config::instance().has("verbose"))
					std::cerr << error << '\n';
				break;
			}
		}

		return result;
	};

	if (url.get_scheme() == "https")
	{
		boost::asio::ssl::context ctx(ssl::context::tls);

		ctx.set_default_verify_paths();
		ctx.set_options(ssl::context::default_workarounds);
		ctx.load_verify_file("/etc/ssl/certs/ca-certificates.crt");

		ssl_socket sock(io_context, ctx);
		boost::asio::connect(sock.lowest_layer(), endpoints);
		sock.lowest_layer().set_option(tcp::no_delay(true));

		(void)SSL_set_tlsext_host_name(sock.native_handle(), host.c_str());

		// Perform SSL handshake and verify the remote host's certificate.
		sock.set_verify_mode(ssl::verify_peer);
#if (BOOST_VERSION / 100 % 1000) >= 73
		sock.set_verify_callback(ssl::host_name_verification(host));
#else
		sock.set_verify_callback(ssl::rfc2818_verification(host));
#endif
		sock.handshake(ssl_socket::client);

		boost::asio::write(sock, req_buffer);

		return reader(sock);
	}
	else
	{
		tcp::socket sock(io_context);
		boost::asio::connect(sock, endpoints);

		boost::asio::write(sock, req_buffer);

		return reader(sock);
	}
}

zeep::http::reply head_request(const zeep::uri &url, std::vector<zeep::http::header> headers)
{
	if (auto &scheme = url.get_scheme(); scheme != "http" and scheme != "https")
		return {};

	// prepare a request
	zeep::http::request req{ "HEAD", url, { 1, 0 }, std::move(headers) };

	return send_request(req, url);
}

zeep::http::reply simple_request(const zeep::uri &url, std::vector<zeep::http::header> headers)
{
	if (auto &scheme = url.get_scheme(); scheme != "http" and scheme != "https")
		return {};

	// prepare a request
	zeep::http::request req{ "GET", url, { 1, 0 }, std::move(headers) };

	return send_request(req, url);
}

zeep::http::reply post_request(const zeep::uri &url, std::vector<zeep::http::header> headers, const std::string &payload)
{
	if (auto &scheme = url.get_scheme(); scheme != "http" and scheme != "https")
		return {};

	// prepare a request

	zeep::http::request req{ "POST", url, { 1, 0 }, std::move(headers) };

	if (auto ct = req.get_header("Content-Type"); ct.empty())
		req.set_content(payload, "application/plain");
	else
	 	req.set_content(payload, ct);

	return send_request(req, url);
}