/**
 * Example file for coela.js, also used for testing
 * 
 */

/*jshint esversion: 6 */
/*jshint node: true */
const TrapSampler = require('./trapdoor.js');
const tsUtils = require("./trapSamplerUtils.js");
const LWESampler = require('./lweSampler.js');
const GGH15Sampler = require('./ggh15Sampler.js');
const matUtils = require('./matUtils.js');
const math = require('mathjs');

class TrapdoorExample {
	constructor(q, n, m, sigma)	{
		this.ts = new TrapSampler(q,n,m,sigma);
	}

	generate() {
		let mats = this.ts.genMatrixWithTrapdoor();
		let U = math.matrix(math.round(math.multiply(math.random([this.ts.n, this.ts.elen]),10)));
		let preImg = this.ts.matPreImage(mats.AA,mats.trap,U);

		return {mats: mats, U: U, preImg: preImg};
	}

	check(res) {
		let AA = res.mats.AA;
		let preImg = res.preImg;
		let U = res.U;
		let rowsU = matUtils.rowSize(U);
		let colsU = matUtils.colSize(U);

		return this.ts.checkMatrixPreImage(AA,preImg,U,rowsU,colsU);
	}
}

class LWEExample {
	constructor(q, n, m, sigma) {
		this.ls = new LWESampler(q, n, m, sigma);
		this.ts = this.ls.ts;
	}

	generate() {
		let mats1 = this.ts.genMatrixWithTrapdoor();
		let mats2 = this.ts.genMatrixWithTrapdoor();
		let lwe = this.ls.shortSecretLWE(mats2.AA);
		let B = lwe.B;
		let preImg = this.ts.matPreImage(mats1.AA,mats1.trap,B);

		return {mats1: mats1, mats2: mats2, B: B, preImg: preImg};
	}

	check(res) {
		let AA0 = res.mats1.AA;
		let preImg = res.preImg;
		let B = res.B;
		let rowsB = matUtils.rowSize(B);
		let colsB = matUtils.colSize(B);

		return this.ts.checkMatrixPreImage(AA0,preImg,B,rowsB,colsB);
	}
}

class GGH15Example {
	constructor(q, n, m, sigma, kappa, branches, bookends) {
		this.ggh15 = new GGH15Sampler(q, n, m, sigma, kappa, branches, bookends);
		this.ls = this.ggh15.lwe;
		this.ts = this.ggh15.ts;
	}

	generate(width) {
		let nodes = this.ggh15.generateNodes();
		let encodings = this.ggh15.generateEdges(nodes.nodes, nodes.bookendMats, width);
		return {nodes: nodes, encodings: encodings};
	}

	check(res,branches,width) {
		let encodings = res.encodings.encs[0];
		let nodes = res.nodes.nodes[0];
		let values = res.encodings.values[0];
		let mult = encodings[0][0];
		let multVals = values[0][0];
		for (let i=1;i<encodings.length;i++) {
			mult = math.multiply(mult,encodings[i][0]);
			multVals = math.multiply(multVals,values[i][0]);
		}

		let AA = nodes[0].AA;
		let AA2 = nodes[2].AA;
		let zt = math.multiply(AA,mult);
		let chk = math.multiply(multVals,AA2);
		let result = tsUtils.matBalance(math.subtract(zt,chk),q);
		let rows = matUtils.rowSize(result);
		let cols = matUtils.colSize(result);
		for (let i=0;i<rows;i++) {
			for (let j=0;j<cols;j++) {
				let entry = result.subset(math.index(i,j));
				if (math.abs(entry) > q/8) {
					return false;
				}
			}
		}
		return true;
	}
}

let q = 1019;
let n = 4;
let m = 3; 
let sigma = 1.0;
let doTrap = false;
let doLWE = false;
let doGGH15 = true;

if (doTrap) {
	let trap = new TrapdoorExample(q, n, m, sigma);
	let resTrap = trap.generate();
	if (!trap.check(resTrap)) {
		throw new Error("Matrix pre-image check for Trapdoor example failed.");
	}
}

if (doLWE) {
	let lwe = new LWEExample(q, n, m, sigma);
	let resLWE = lwe.generate();
	if (!lwe.check(resLWE)) {
		throw new Error("Matrix pre-image check for LWE example failed.");
	}
}

if (doGGH15) {
	let kappa = 3;
	// for (let branches=1;branches<3;branches++) {
	// 	for (let width=1;width<3;width++) {
	// 		let ggh15NoBookends = new GGH15Example(q, n, m, sigma, kappa, branches, false);
	// 		let resGGH15NoBookends = ggh15NoBookends.generate(width);

	// 		let ggh15Bookends = new GGH15Example(q, n, m, sigma, kappa, branches, true);
	// 		let resGGH15Bookends = ggh15Bookends.generate(width);
	// 	}
	// }

	let ggh15 = new GGH15Example(q, n, m, sigma, kappa, 1, false);
	let resGGH15 = ggh15.generate(1);
	if (!ggh15.check(resGGH15,1,1)) {
		throw new Error("GGH15 check failed");
	}
}



