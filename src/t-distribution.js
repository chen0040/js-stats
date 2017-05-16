var jsstats = jsstats || {};

(function(jss){
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

	jsstats.TDistribution = TDistribution;
})(jsstats);

if(module){
	module.exports = jsstats.TDistribution;
}