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
		var Z = this.Sqrt2 * this.invErrorFunc(2 * p - 1);
        return Z * this.sd + this.mean;
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

	jss.TDistribution = TDistribution;

    var FDistribution = function(df1, df2) {
        this.df1 = df1;
        this.df2 = df2;
        this.EPSILON = 0.0000000001;
    };

    FDistribution.prototype.L504 = function(a, f, b, iv)
    {
        var q = a * f / (a * f + b);
        var sa = Math.sqrt(q);
        var sl = Math.log(sa);
        var ca = Math.sqrt(1 - q);
        var cl = Math.log(ca);
        var al = Math.atan(sa / Math.sqrt(-sa * sa + 1));
        var fp = 1 - 2 * al / Math.PI;
        var r = 0.0;
        if (b != 1)
        {
            var c = Math.log(2 * sa / Math.PI);
            fp -= Math.exp(c + cl);
            if (b != 3)
            {
                var n = Math.floor((b - 3) / 2);
                for (var i = 1; i <= n; i++)
                {
                    var x = 2 * i + 1;
                    r += Math.log((x - 1) / x);
                    var rr = r + cl * x + c;
                    if (rr > -78.4)
                    {
                        fp -= Math.exp(rr);
                    }
                }
            }
        }

        if (a != 1)
        {
            var c = r;

            if (b > 1)
            {
                c += Math.log(b - 1);
            }

            c += Math.log(2 / Math.PI) + sl + cl * b;

            if (c > -78.4) { fp += Math.exp(c); }

            if (a != 3)
            {
                var n = Math.floor((a - 3) / 2);
                r = 0;
                for (var i = 1; i <= n; i++)
                {
                    var x = i * 2 + 1;
                    r += Math.log((b + x - 2) / x);
                    var rr = r + sl * (x - 1) + c;
                    if (rr > -78.4) { fp += Math.exp(rr); }
                }
            }
        }
        return fp;

    };

    FDistribution.prototype.L401 = function(a, f, b, iv)
    {
        var q = a * f / (a * f + b);
        var ql = Math.log(q);
        var fp = 0.0;
        var c = Math.log(1 - q) * b / 2;
        if (c > -78.4)
        {
            fp = Math.exp(c);
        }

        if (a != 2)
        {
            var n = Math.floor(a / 2 - 1);
            var r = 0.0;
            for (var i = 1; i <= n; i++)
            {
                var x = 2 * i;
                r += Math.log(b + x - 2) - Math.log(x) + ql;
                if (r + c > -78.4)
                {
                    fp += Math.exp(r + c);
                }
            }
        }

        if (iv == 1)
        {
            fp = 1 - fp;
        }

        return fp;
    };

    FDistribution.prototype.ProbF = function(dn, dd, fr)
    {
        var f = fr;
        var a = dn;
        var b = dd;
        var iv = 0;

        if (Math.floor(a / 2) * 2 == a)
        {
            //even numerator df
            var fp = this.L401(a, f, b, iv);
            return fp;
        }
        else if (Math.floor(b / 2) * 2 != b)
        {
            var fp = this.L504(a, f, b, iv);
            return fp;
        }

        f = 1 / f;
        a = dd;
        b = dn;
        iv = 1;
        return this.L401(a, f, b, iv);

    };

    FDistribution.prototype.cumulativeProbability = function(F) {
        if (this.df1 > .01 & this.df2 > .01 & F > this.EPSILON)
        {
            var p = 1 - this.ProbF(this.df1, this.df2, F);
            return p;
        }
        else
        {
            console.error("df1, df2, and F must be numbers greater than 0.");
        }
    };

    jss.FDistribution = FDistribution;

    var ChiSquareDistribution = function(df) {
        this.df = df;
    };

    ChiSquareDistribution.prototype.ChiSquaredProbability = function(x) {
        var a, y = 0, s, e, c, z, val;
        var df = this.df;
        var bigx = 20.0;
        var logSqrtPi = Math.log(Math.sqrt(Math.PI));
        var rezSqrtPi = 1 / Math.sqrt(Math.PI);
        if (x <= 0 || df < 1)
            return (1);
        a = 0.5 * x;
        even = ((parseInt(2 * (df / 2), 2)) == df);
        if (df > 1)
            y = Math.exp(-a); //((-a < -bigx) ? 0.0 : Math.exp (-a));
        s = (even ? y : (2.0 * (new jsstats.NormalDistribution(0.0, 1.0).cumulativeProbability(-Math.sqrt(x)))));
        if (df > 2)
        {
            x = 0.5 * (df - 1.0);
            z = (even ? 1.0 : 0.5);
            if (a > bigx)
            {
                e = (even ? 0.0 : logSqrtPi);
                c = Math.log(a);
                while (z <= x)
                {
                    e = Math.log(z) + e;
                    val = c * z - a - e;
                    s += Math.exp(val); //((val < -bigx) ? 0.0 : Math.exp (val));
                    z += 1.0;
                }
                return (s);
            }
            else
            {
                e = (even ? 1.0 : (rezSqrtPi / Math.sqrt(a)));
                c = 0.0;
                while (z <= x)
                {
                    e = e * (a / z);
                    c = c + e;
                    z += 1.0;
                }
                return (c * y + s);
            }
        }
        else
        {
            return (s);
        }
    };

    ChiSquareDistribution.prototype.cumulativeProbability = function(x) {
        return 1 - this.ChiSquaredProbability(x);
    };

    jss.ChiSquareDistribution = ChiSquareDistribution;

})(jsstats);

var module = module || {};
if(module) {
	module.exports = jsstats;
}