var jsstats = jsstats || {};

(function(jss){
	jss.NormalDistribution = function(mean, sd){
		if(!mean) {
			mean = 0.0;
		}
		if(!sd) {
			sd = 1.0;
		}
		this.mean = mean;
		this.sd = sd;
	};
})(jsstats);

if(module) {
	module.exports = jsstats.NormalDistribution;
}