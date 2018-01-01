/**
 * Sampler for sampling pre-images with respect to LWE samples
 * 
 */

/*jshint esversion: 6 */
/*jshint node: true */
'use strict';
const TrapSampler = require('./trapdoor-sampler.js');
const matUtils = require('./mat-utils.js');
const math = require('mathjs');

class LWESampler {
	constructor(q, n, m, sigma) {
		this.ts = new TrapSampler(q, n, m, sigma);
	}

	/**
	 * Samples an LWE instance with a short secret
	 * @param {Matrix | Array} A [Uniform matrix]
	 * @return {Matrix | Array literal} [LWE sample]
	 */
	shortSecretLWE(A) {
		let ts = this.ts;
		let S = matUtils.genShortMatrix(ts.n,ts.n,ts.dist);
		let B = this.getLWEMatrix(A,S);
		return {A: A, B: B};
	}

	/**
	 * Samples an LWE sample with respect to a specific secret
	 * @param {Matrix | Array} A 
	 * @param {Matrix | Array} S 
	 * @returns {Matrix | Array literal}
	 */
	specificSecretLWE(A,S) {
		let ts = this.ts;
		let rowsS = matUtils.rowSize(S);
		let colsS = matUtils.colSize(S);
		if (rowsS != ts.n || colsS != ts.n) {
			throw new Error("Specified matrix S: " + S + " does not have correct dimensions, rows: " + rowsS + "; cols: " + colsS + ".");
		}
		let B = this.getLWEMatrix(A,S);
		return {A: A, B: B};
	}

	/**
	 * Computes AS+E for a randomly sampled error E
	 * @param {Matrix | Array} A 
	 * @param {Matrix | Array} S
	 * @returns {Matrix | Array literal}
	 */
	getLWEMatrix(A,S) {
		let ts = this.ts;
		let E = matUtils.genShortMatrix(ts.n, ts.elen, ts.dist);
		let SA = math.multiply(S,A);
		let B = math.mod(math.add(SA, E), ts.q);
		return B;
	}
}

module.exports = LWESampler;
