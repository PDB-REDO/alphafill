const { CleanWebpackPlugin } = require('clean-webpack-plugin');
const path = require('path');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');

const SCRIPTS = path.resolve(__dirname, "webapp");
const SCSS = __dirname + "/scss/";
const DEST = path.resolve(__dirname, "docroot");

module.exports = (env) => {

	const PRODUCTION = env != null && env.PRODUCTION;

	const webpackConf = {

		entry: {
			'pdb-redo-bootstrap': path.resolve(SCSS, "pdb-redo-bootstrap.scss"),
			'index': path.resolve(SCRIPTS, "index.js"),
			'model': path.resolve(SCRIPTS, "model.js"),
			'optimized': path.resolve(SCRIPTS, "optimized.js"),
			'lists': path.resolve(SCRIPTS, "lists.js"),
			'wait': path.resolve(SCRIPTS, "wait.js"),
		},

		output: {
			path: DEST,
			crossOriginLoading: 'anonymous',
			filename: "scripts/[name].js"
		},

		module: {
			rules: [
				{
					test: /\.js/,
					exclude: /node_modules/,
					use: {
						loader: "babel-loader",
						options: {
							presets: ['@babel/preset-env']
						}
					}
				},

				{
					test: /\.(sa|sc|c)ss$/i,
					use: [
						MiniCssExtractPlugin.loader,
						"css-loader",
						"postcss-loader",
						"sass-loader"
					]
				},

				{
					test: /\.woff(2)?(\?v=[0-9]\.[0-9]\.[0-9])?$/,
					include: path.resolve(__dirname, './node_modules/bootstrap-icons/font/fonts'),
					type: 'asset/resource',
					generator: {
						filename: 'fonts/[name][ext]'
					}
				}
			]
		},

		resolve: {
			extensions: ['.js', '.scss'],
		},

		optimization: { minimizer: [] },

		target: 'web',

		plugins: [
			new MiniCssExtractPlugin({
				filename: "css/[name].css"
			})
		]
	};

	if (PRODUCTION) {
		webpackConf.mode = "production";

		webpackConf.plugins.push(
			new CleanWebpackPlugin({
				verbose: true,
				cleanOnceBeforeBuildPatterns: [
					'css/*',
					'fonts/*',
					'scripts/*',
				]
				
			})
		);
	} else {
		webpackConf.mode = "development";
		webpackConf.devtool = 'source-map';
	}

	return webpackConf;
};
