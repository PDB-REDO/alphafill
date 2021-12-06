const MiniCssExtractPlugin = require('mini-css-extract-plugin');
const { CleanWebpackPlugin } = require('clean-webpack-plugin');
const webpack = require('webpack');

const SCRIPTS = __dirname + "/webapp/";
const SCSS = __dirname + "/scss/";
const DEST = __dirname + "/docroot/";

module.exports = {
	mode: "development",
	entry: {
		'pdb-redo-bootstrap': SCSS + "pdb-redo-bootstrap.scss",
		'molstar': SCRIPTS + "molstar.tsx",
		'index': SCRIPTS + "index.js",
		'model': SCRIPTS + "model.js",
		'lists': SCRIPTS + "lists.js",
	},

	output: {
		path: DEST,
		filename: "./scripts/[name].js"
	},

	devtool: "source-map",

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
				test: /\.tsx?$/,
				exclude: /node_modules/,
				use: {
					loader: "babel-loader",
					options: {
						presets: ['@babel/preset-env']
					}
				}
			},
			{
				test: /\.css$/,
				use: [
					// 'style-loader',
					MiniCssExtractPlugin.loader,
					'css-loader'
				]
			},
			{
				test: /\.(eot|svg|ttf|woff(2)?)(\?v=\d+\.\d+\.\d+)?/,
				loader: 'file-loader',
				options: {
					name: '[name].[ext]',
					outputPath: 'fonts/',
					publicPath: '../fonts/'
				}
			},
			{
				test: /\.s[ac]ss$/i,
				use: [
					MiniCssExtractPlugin.loader,
					'css-loader',
					'sass-loader'
				]
			},
			{
				test: /\.(png|jpg|gif)$/,
				use: [
					{
						loader: 'file-loader',
						options: {
							outputPath: "css/images",
							publicPath: "images/"
						},
					},
				]
			}
		]
	},


	resolve: {
		extensions: ['.tsx', '.ts', '.js'],
	},

	plugins: [
		new CleanWebpackPlugin({
			cleanOnceBeforeBuildPatterns: [
				'css/**/*',
				'css/*',
				'scripts/**/*',
				'fonts/**/*'
			]
		}),
		new webpack.ProvidePlugin({
			$: 'jquery',
			jQuery: 'jquery'
		}),
		new MiniCssExtractPlugin({
			filename: './css/[name].css',
			chunkFilename: './css/[id].css'
		})
	]
};
