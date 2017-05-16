var expect    = require("chai").expect;
var NormalDistribution = require("../src/normal-distribution");

describe("Create normal distribution", function() {
  describe("default constructor", function() {
    it("has mean of 0.0 and sd of 1.0", function() {
    	var distribution = new NormalDistribution();
    	expect(distribution.mean).to.equal(0.0);
    	expect(distribution.sd).to.equal(1.0);
    });
  });

  describe("Constructor with arguments", function() {
    it("has user-defined mean and sd", function() {
    	var distribution = new NormalDistribution(5.0, 12.0);
    	expect(distribution.mean).to.equal(5.0);
    	expect(distribution.sd).to.equal(12.0);
    });
  });
});