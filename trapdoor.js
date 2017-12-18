/**
 * Trapdoor javascript file for instantiating different sampling methods
 * 
 */

/*jshint esversion: 6 */
/*jshint node: true */
'use strict';
const rye = require('rye');
const math = require('mathjs');
const gaussian = require('gaussian');
const sjcl = require('sjcl');
const tsUtils = require('./trapSamplerUtils.js');
const matUtils = require('./matUtils.js');
const ringUtils = require('./ringUtils.js');
// const Polynomial = require('./poly.js');
const PrimeField = rye.PrimeField;

class TrapSampler {
    constructor(q, n, m, sigma) {
        this.q = q;
        this.n = n;
        this.m = m;
        this.l = math.ceil(math.log(q, 2));
        this.ZZq = new PrimeField(q);
        this.dist = new DiscreteGaussian(0, sigma);
        this.g = tsUtils.getGadgetVector(this.l);
        this.G = tsUtils.getGadget(this.g, this.n);
        this.zMap = tsUtils.createBuckets(this.q, this.l, this.g, this.dist, null);
        this.elen = this.m + this.n * this.l;
        this.trapWidth = this.n * this.l;
    }

    /**
     * Generates a uniform matrix AA along with a trapdoor matrix RR
     * @return {Array} [AA,RR]
     */
    genMatrixWithTrapdoor() {
        let A = math.mod(math.randomInt([this.n, this.m], this.q), this.q);
        let R = matUtils.genShortMatrix(this.m, this.trapWidth, this.dist);
        let A1 = tsUtils.balance(math.subtract(this.G, math.multiply(A, R)), this.q);
        let IR = math.eye(this.trapWidth);
        // AA \in \ZZ_q^{n x (elen)}
        let AA = matUtils.augment(A, A1);
        // RR \in \ZZ_q^{(elen) x trapWidth}
        let RR = matUtils.stack(R, IR);
        return [AA, RR];
    }

    /**
     * * Return a pre-image matrix X such that AA*X = U mod q for a vector U
     * @param {Matrix | Array} AA 
     * @param {Matrix | Array} RR 
     * @param {Matrix | Array} U 
     */
    matPreImage(AA, RR, U) {
        let X = [];
        for (let i = 0; i < matUtils.rowSize(U); i++) {
            let x = this.vecPreImage(AA, RR, matUtils.getRow(U, i));
            X.push(x);
        }
        let matX = math.transpose(tsUtils.matBalance(math.matrix(X), this.q));
        if (!math.equal(math.mod(math.multiply(AA, matX), this.q), U)) {
            throw new Error("Matrices do not equal after pre-image computation AA*X:" + math.mod(math.multiply(AA, matX), this.q).valueOf() + ", U:" + U.valueOf());
        }
        return matX;
    }

    /**
     * Return a pre-image x such that AA*x = u mod q for a vector u
     * @param {Matrix | Array} AA 
     * @param {Matrix | Array} RR 
     * @param {Vector} u (n x 1)
     * @return {Vector} x
     */
    vecPreImage(AA, RR, u) {
        let pArr = [];
        for (let i = 0; i < this.elen; i++) {
            pArr.push(this.dist.sample);
        }
        let p = matUtils.wrapVector(pArr);
        // v = (u-AA*p) mod q
        let v = (math.mod(math.subtract(u, math.multiply(AA, p)), this.q));
        let z = null;
        try {
            z = this.zPreImg(v);
        } catch (error) {
            console.error(error);
            return null;
        }
        // x = (p + RR*z) mod q 
        let x = math.mod(math.add(p, math.multiply(RR, z)), this.q);
        // res = (AA*x) mod q
        let res = math.mod(math.multiply(AA, x), this.q);
        if (!math.equal(res, u)) {
            throw new Error("Incorrect vector pre-image calculated, res: " + res + ", u: " + u);
        } else {
            return x;
        }
    }

    /**
     * For a given vector v, constructs a z such that 
     * G*z = v mod q
     * @param {Vector (mod q)} v
     * @return {Vector} z 
     */
    zPreImg(v) {
        let z = [];
        let vArr = v.valueOf();
        let self = this;
        vArr.forEach(function (element) {
            let samps = self.zMap[element];
            if (matUtils.colSize(samps) != 0) {
                z = z.concat(samps.pop().valueOf());
            } else {
                console.log("There are no samples left for: ", + element);
                self.zMap = tsUtils.createBuckets(self.q, self.l, self.g, self.dist, self.zMap);
                throw new Error("New z values sampled, please retry.");
            }
        });
        return matUtils.wrapVector(z);
    }
}

class DiscreteGaussian {
    constructor(mean, variance) {
        this.dist = gaussian(mean, variance);
    }
    get sample() {
        return math.round(this.dist.ppf(math.random()));
    }
}