/**
 * Trapdoor sampler utils
 */
const math = require('mathjs');
const matUtils = require('./matUtils.js');

const methods = {
    matBalance: function(X,q) {
        let newX = [];
        X.forEach(function (value, index) {
            let rowIndex = math.subset(index, math.index(0));
            let colIndex = math.subset(index, math.index(1));
            X.subset(math.index(rowIndex, colIndex), methods.balance(value,q));
        });
        return X;
    },

    /**
     * Writes the value x mod q into the 
     * the bounds -q/2 < x < q/2
     * @param {Integer} x 
     * @param {Prime number} q 
     * @return {Number} xq
     */
    balance: function(x,q) {
        if (!q) {
            throw new Error("bad error");
        }
        let xq = math.mod(x,q);
        if (xq > math.floor(q/2)) {
            xq = xq - q;
        }
        return xq
    },

    /**
     * Returns an l-dimensional gadget vector
     * @param {Dimension} l 
     * @return {Vector} g
     */
    getGadgetVector: function(l) {
        let arr = [];
        for (let i=0; i<l; i++) {
            arr.push(math.pow(2,i));
        }
        return matUtils.wrapVector(arr)
    },

    /**
     * Returns an (n x ln) gadget matrix
     * @param {Gadget Vector} g 
     * @param {Dimension} n 
     * @return {Matrix} G
     */
    getGadget: function(g,n) {
        return matUtils.tensorProduct(math.eye(n), g);
    },

    /**
     * Creates a set of buckets for sampling z values
     * @param {Prime} q 
     * @param {Gadget Vector dimension} l 
     * @param {Gadget Vector} g 
     * @param {Discrete Gaussian} dist 
     * @return {Buckets of DG samples} zMap
     */
    createBuckets: function(q,l,g,dist,zMap) {
        if (!zMap) {
            zMap = new Map();
        }
        for (let i=0; i<q*l; i++) {
            zVals = [];
            for (let j=0; j<l; j++) {
                zVals.push(dist.sample);
            }
            z = matUtils.wrapVector(zVals);
            zq = math.mod(math.multiply(g,z),q);
            if (!zMap[zq]) {
                zMap[zq] = [];
            }
            zMap[zq] = zMap[zq].concat(z);
        }
        return zMap;
    },
};

module.exports = methods;