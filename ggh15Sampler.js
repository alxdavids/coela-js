/**
 * Sampler for sampling pre-images with respect to LWE samples
 * 
 */

/*jshint esversion: 6 */
/*jshint node: true */
'use strict';
const TrapSampler = require('./trapdoor.js');
const LWESampler = require('./lweSampler.js');
const matUtils = require('./matUtils.js');
const math = require('mathjs');

class GGH15Sampler {
	constructor(q, n, m, sigma, kappa, branches, bookends) {
		this.lwe = new LWESampler(q, n, m, sigma);
		this.ts = this.lwe.ts;
		// 0 < n = number of branches in the graph
		this.branches = branches;
		// Number of nodes in each branch
		this.kappa = kappa;
		// Specify whether bookend matrices are required for the graph
		this.bookends = bookends;
	}

	/**
	 * Generates the nodes for a graph
	 * @returns {Matrix | Array literal} [Graph nodes]
	 */
	generateNodes() {
		if (this.branches < 1) {
			throw new Error("Graph must have at least one branch");
		}

		let kappa = this.kappa;
		let ts = this.ts;
		// Generate trapdoor matrices for each node in the graph
		let nodes = [];
		for (let i=0;i<this.branches;i++) {
			nodes[i] = [];
			for (let j=0;j<kappa;j++) {
				let mats = ts.genMatrixWithTrapdoor();
				nodes[i][j] = mats;
			}
		}

		// Also sample bookend matrices, if need be
		let bookendMats;
		if (this.bookends) {
			let start = ts.genMatrixWithTrapdoor();
			let end = ts.genMatrixWithTrapdoor();
			bookendMats = {start: start, end: end};
		}
		return {bookendMats: bookendMats, nodes: nodes};
	}

	/**
	 * Generate encodings for graph edges
	 * @param {Matrix | Array} nodes 
	 * @param {Matrix | Array} bookendMats 
	 * @param {Matrix | Array} width
	 * @returns {GGH15 literal}
	 */
	generateEdges(nodes,bookendMats,width) {
		if (width < 1) {
			throw new Error("Width of GGH15 sampler must be greater than 0");
		}
		let ts = this.ts;
		let lwe = this.lwe;
		let branches = this.branches;
		let kappa = this.kappa;
		let encodings = [];
		let values = [];
		let bookendEncodings = [];
    	let bookendVals = [];
		for (let i=0;i<branches;i++) {
			encodings[i] = [];
			values[i] = [];
			for (let j=0;j<kappa-1;j++) {
				let matsSource = nodes[i][j];
				let matsSink = nodes[i][j+1];
                let res = this.computePaths(matsSource, matsSink, width);
				encodings[i][j] = res.encs;
				values[i][j] = res.vals;
			}

			if (this.bookends) {
				let matsSourceStart = bookendMats.start;
				let matsSinkStart = nodes[i][0];
				let matsSourceEnd = nodes[i][nodes[i].length-1];
				let matsSinkEnd = bookendMats.end;
				let resStart = this.computePaths(matsSourceStart, matsSinkStart, width);
				let resEnd = this.computePaths(matsSourceEnd,matsSinkEnd,width);
				bookendEncodings[i] = {start: resStart.encs, end: resEnd.encs};
        		bookendVals[i] = {start: resStart.vals, end: resEnd.vals};
			}
		}
		return {encs: encodings, values: values, bookendEncs: bookendEncodings, bookendVals: bookendVals};
	}

	/**
	 * Returns an encoding of S on the path AA0 -> AA1
	 * @param {Matrix | Array} AA0 
	 * @param {Matrix | Array} RR0 
	 * @param {Matrix | Array} AA1 
	 * @param {Matrix | Array} S 
	 */
	generateEncoding(AA0,RR0,AA1,S) {
		let ts = this.ts;
		let lwe = this.lwe;
		let lweSample = lwe.specificSecretLWE(AA1,S);
		let encS = ts.matPreImage(AA0,RR0,lweSample.B);
		return encS;
	}

	/**
	 * Returns a number of encodings corresponding to the designated path
	 * @param {Matrix | Array} matsSource 
	 * @param {Matrix | Array} matsSink 
	 * @returns {Matrix | Array literal}
	 */
	computePaths(matsSource, matsSink, width) {
		let ts = this.ts;
		let ls = this.lwe;
		let encs = [];
		let vals = [];
		for (let w = 0; w < width; w++) {
			let S = matUtils.genShortMatrix(ts.n, ts.n, ts.dist);
			let encS = this.generateEncoding(matsSource.AA, matsSource.trap, matsSink.AA, S);
			encs[w] = encS;
			vals[w] = S;
		}
		return {encs: encs, vals: vals};
	}
}

module.exports = GGH15Sampler;