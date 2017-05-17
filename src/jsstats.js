var jsstats = jsstats || {};

(function(jss){
	var NormalDistribution = function(mean, sd){
		if(!mean) {
			mean = 0.0;
		}
		if(!sd) {
			sd = 1.0;
		}
		this.mean = mean;
		this.sd = sd;
		this.Sqrt2 = 1.4142135623730950488016887;
		this.Sqrt2PI = 2.50662827463100050242E0;
		this.lnconstant = -Math.log(this.Sqrt2PI * sd);
	};

	NormalDistribution.prototype.sample = function() {

	};

	NormalDistribution.prototype.cumulativeProbability = function(x) {
		var z = (x - this.mean) / (this.Sqrt2 * this.sd);
		return 0.5 + 0.5 * this.errorFunc(z);
	};

	NormalDistribution.prototype.invCumulativeProbability = function(p) {
		return this.Sqrt2 * this.invErrorFunc(2 * p - 1);
	};

	NormalDistribution.prototype.errorFunc = function(z){
		var t = 1.0 / (1.0 + 0.5 * Math.abs(z));

	    // use Horner's method
	    var ans = 1 - t * Math.exp(-z * z - 1.26551223 +
	                                        t * (1.00002368 +
	                                        t * (0.37409196 +
	                                        t * (0.09678418 +
	                                        t * (-0.18628806 +
	                                        t * (0.27886807 +
	                                        t * (-1.13520398 +
	                                        t * (1.48851587 +
	                                        t * (-0.82215223 +
	                                        t * (0.17087277))))))))));
            if (z >= 0) return ans;
            else return -ans;
	};

	NormalDistribution.prototype.invErrorFunc = function(x){
		var z;
        var a = 0.147;
        var the_sign_of_x;
        if (0 == x)
        {
            the_sign_of_x = 0;
        }
        else if (x > 0)
        {
            the_sign_of_x = 1;
        }
        else
        {
            the_sign_of_x = -1;
        }

        if (0 != x)
        {
            var ln_1minus_x_sqrd = Math.log(1 - x * x);
            var ln_1minusxx_by_a = ln_1minus_x_sqrd / a;
            var ln_1minusxx_by_2 = ln_1minus_x_sqrd / 2;
            var ln_etc_by2_plus2 = ln_1minusxx_by_2 + (2 / (Math.PI * a));
            var first_sqrt = Math.sqrt((ln_etc_by2_plus2 * ln_etc_by2_plus2) - ln_1minusxx_by_a);
            var second_sqrt = Math.sqrt(first_sqrt - ln_etc_by2_plus2);
            z = second_sqrt * the_sign_of_x;
        }
        else
        { // x is zero
            z = 0;
        }
        return z;
	};

	jss.NormalDistribution = NormalDistribution;

	var TDistribution = function(df){
		if(df){
			this.df = df;
		}
	};

	TDistribution.prototype.LogGamma = function(Z) {
		with (Math) {
			var S=1+76.18009173/Z-86.50532033/(Z+1)+24.01409822/(Z+2)-1.231739516/(Z+3)+.00120858003/(Z+4)-.00000536382/(Z+5);
			var LG= (Z-.5)*log(Z+4.5)-(Z+4.5)+log(S*2.50662827465);
		}
		return LG
	}

	TDistribution.prototype.Betinc = function(X,A,B) {
		var A0=0;
		var B0=1;
		var A1=1;
		var B1=1;
		var M9=0;
		var A2=0;
		var C9;
		while (Math.abs((A1-A2)/A1)>.00001) {
			A2=A1;
			C9=-(A+M9)*(A+B+M9)*X/(A+2*M9)/(A+2*M9+1);
			A0=A1+C9*A0;
			B0=B1+C9*B0;
			M9=M9+1;
			C9=M9*(B-M9)*X/(A+2*M9-1)/(A+2*M9);
			A1=A0+C9*A1;
			B1=B0+C9*B1;
			A0=A0/B1;
			B0=B0/B1;
			A1=A1/B1;
			B1=1;
		}
		return A1/A
	}

	TDistribution.prototype.cumulativeProbability = function(X, df) {
		if(!df) {
			df = this.df;
		}

	    with (Math) {
			if (df<=0) {
				console.error("Degrees of freedom must be positive");
			} else {
				A=df/2;
				S=A+.5;
				Z=df/(df+X*X);
				BT=exp(this.LogGamma(S)-this.LogGamma(.5)-this.LogGamma(A)+A*log(Z)+.5*log(1-Z));
				if (Z<(A+1)/(S+2)) {
					betacdf=BT*this.Betinc(Z,A,.5)
				} else {
					betacdf=1-BT*this.Betinc(1-Z,.5,A)
				}
				if (X<0) {
					tcdf=betacdf/2
				} else {
					tcdf=1-betacdf/2
				}
			}
			tcdf=round(tcdf*100000)/100000;
		}
    	return tcdf;
	};

	TDistribution.prototype.invCumulativeProbability = function(p, df) {
		if(!df){
			df = this.df;
		}
		var delta = 0.005;
        
        if(p >= 0.5) {
            var Z1 = 0;
            for(Z = 0; Z < 100; Z++) {
                if(this.cumulativeProbability(Z, df) >= p){
                    break;
                }
                Z1 = Z;
            }
            var Z2 = Z1;
            for(var Z = 0.0; Z < 100.0; Z+=1.0) {
                
                if(this.cumulativeProbability(Z1 + Z / 100.0) >= p){
                    break;
                }
                Z2 = Z1 + (Z)/100.0;
            }
            var Z3 = Z2;
            for(var Z = 0.0; Z < 100.0; Z+=1.0) {
                
                if(this.cumulativeProbability(Z2 + Z / 10000.0) >= p){
                    break;
                }
                Z3 = Z2 + (Z)/10000.0;
            }
            return Z3;
        } else {
            var Z1 = 0;
            for(var Z = 0; Z < 100; Z++) {
                if(this.cumulativeProbability(-Z, df) <= p){
                    break;
                }
                Z1 = Z;
            }
            var Z2 = Z1;
            for(var Z = 0.0; Z < 100.0; Z+=1.0) {
                
                if(this.cumulativeProbability(-Z1 - Z / 100.0) <= p){
                    break;
                }
                Z2 = Z1 + (Z) / 100.0;
            }
            var Z3 = Z2;
            for(var Z = 0.0; Z < 100.0; Z+=1.0) {
                
                if(this.cumulativeProbability(-Z2 - Z / 10000.0) <= p){
                    break;
                }
                Z3 = Z2 + (Z)/10000.0;
            }
            return -Z3;
        }
    };

	jsstats.TDistribution = TDistribution;


})(jsstats);

if(module) {
	module.exports = jsstats;
}