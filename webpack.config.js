const { CleanWebpackPlugin } = require('clean-webpack-plugin');
const path = require('path');
const MiniCssExtractPlugin = require('mini-css-extract-plugin');

const SCRIPTS = __dirname + "/webapp/";
const SCSS = __dirname + "/scss/";
const DEST = __dirname + "/docroot/scripts/";

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
			crossOriginLoading: 'anonymous'
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
						/* PRODUCTION ?  */MiniCssExtractPlugin.loader/*  : "style-loader" */,
						"css-loader",
						"postcss-loader",
						"sass-loader"
					]
				},

				{
					test: /\.woff(2)?(\?v=[0-9]\.[0-9]\.[0-9])?$/,
					include: path.resolve(__dirname, './node_modules/bootstrap-icons/font/fonts'),
					type: 'asset/resource'
				}
			]
		},


		resolve: {
			extensions: ['.js', '.scss'],
		},

		plugins: [
			new MiniCssExtractPlugin({}),
			new CleanWebpackPlugin({
				verbose: true,
				cleanOnceBeforeBuildPatterns: [
					'css/**/*',
					'css/*',
					'scripts/**/*',
					'fonts/**/*'
				]
			})
		],

		optimization: {
			minimizer: []
		}
	};

	if (PRODUCTION) {
		webpackConf.mode = "production";

		// webpackConf.plugins.push(
		// 	new CleanWebpackPlugin({
		// 		verbose: true
		// 	})/* ,
		// 	new MiniCssExtractPlugin({}) */
		// );
	} else {
		webpackConf.mode = "development";
		webpackConf.devtool = 'source-map';
	}

	return webpackConf;
};

