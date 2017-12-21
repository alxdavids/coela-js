/**
 * Sampler for sampling pre-images with respect to LWE samples
 * 
 */

/*jshint esversion: 6 */
/*jshint node: true */
'use strict';
const TrapSampler = require('./trapdoor.js');
const matUtils = require('./matUtils.js');
const math = require('mathjs');

class LWESampler {
	constructor(q, n, m, sigma) {
		this.ts = new TrapSampler(q, n, m, sigma);
	}

	/**
	 * Samples an LWE instance with a short secret
	 * @param  {Matrix | Array} A [Uniform matrix]
	 * @return {Matrix | Array literal} [LWE sample]
	 */
	sampleShortLWE(A) {
		let ts = this.ts;
		let S = matUtils.genShortMatrix(ts.elen,ts.n,ts.dist);
		let E = matUtils.genShortMatrix(ts.n,ts.n,ts.dist);
		let AS = math.multiply(A,S);
		let B = math.mod(math.add(AS,E), ts.q);
		return {A: A, B: B};
	}
}

module.exports = LWESampler;
