/**
 * Ring utility functions
 */

/*jshint esversion: 6 */
/*jshint node: true */
 'use strict';
const rye = require("rye");
const PolynomRing = rye.PolynomRing;
const FactorRing = rye.FactorRing;

const methods = {
    initFactorRing: function(field, dim) {
        let PR = new PolynomRing(field);
        let quotient = methods.getQuotient(dim);
        let FR = new FactorRing(PR, PR.polynom(quotient));
        return FR;
    },

    getQuotient: function(dim) {
        let quo = [];
        for (let i=0; i<dim+1; i++) {
            if (i == 0 || i == dim) {
                quo.push(1);
            } else {
                quo.push(0);
            }
        }
        return quo;
    }
}
module.exports = methods;