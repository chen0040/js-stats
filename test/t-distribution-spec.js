var expect    = require("chai").expect;
var jsstats = require("../src/jsstats");

describe("Create t distribution", function() {
  describe("default constructor", function() {
    var distribution = new jsstats.TDistribution(10.0);
    it("has df of 10.0", function() {
    	expect(distribution.df).to.equal(10.0); 
      
    });
  });

  describe('run cumulative probability', function(){
    it('has probability of 0.5 at t_df = 0', function(){
      var distribution = new jsstats.TDistribution(10.0);
      expect(distribution.cumulativeProbability(0.0)).to.equal(0.5);
    });
    it('has t_df = 0 with probability of 0.5', function(){
      var distribution = new jsstats.TDistribution(10);
      expect(distribution.invCumulativeProbability(0.5)).to.equal(0.0);
    });
  });


});