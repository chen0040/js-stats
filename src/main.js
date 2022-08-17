

export class NormalDistribution { 

    constructor(mean, sd) {
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
	}

	sample() {
	}

	cumulativeProbability(x) {
		const z = (x - this.mean) / (this.Sqrt2 * this.sd);
		return 0.5 + 0.5 * this.errorFunc(z);
	}

	invCumulativeProbability(p) {
		const Z = this.Sqrt2 * this.invErrorFunc(2 * p - 1);
        return Z * this.sd + this.mean;
	}

	errorFunc(z) {
		const t = 1.0 / (1.0 + 0.5 * Math.abs(z));

	    // use Horner's method
	    const ans = 1 - t * Math.exp(-z * z - 1.26551223 +
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
	}

	invErrorFunc(x) {
		let z;
        const a = 0.147;
        let the_sign_of_x;
        if (0 == x) {
            the_sign_of_x = 0;
        } else if (x > 0) {
            the_sign_of_x = 1;
        } else {
            the_sign_of_x = -1;
        }

        if (0 != x) {
            const ln_1minus_x_sqrd = Math.log(1 - x * x);
            const ln_1minusxx_by_a = ln_1minus_x_sqrd / a;
            const ln_1minusxx_by_2 = ln_1minus_x_sqrd / 2;
            const ln_etc_by2_plus2 = ln_1minusxx_by_2 + (2 / (Math.PI * a));
            const first_sqrt = Math.sqrt((ln_etc_by2_plus2 * ln_etc_by2_plus2) - ln_1minusxx_by_a);
            const second_sqrt = Math.sqrt(first_sqrt - ln_etc_by2_plus2);
            z = second_sqrt * the_sign_of_x;
        } else { // x is zero
            z = 0;
        }
        return z;
	}
}

export class TDistribution {

    constructor(df) {
		if (df){
			this.df = df;
		}
	}

	LogGamma(Z) {		
        const S = 1 + 76.18009173/Z - 86.50532033/(Z + 1) + 24.01409822/(Z + 2) 
            - 1.231739516/(Z + 3) + .00120858003/(Z + 4) - .00000536382/(Z + 5);
        const LG = (Z - .5)*Math.log(Z + 4.5) - (Z + 4.5) + Math.log(S*2.50662827465);
		return LG;
	}

	Betinc(X, A, B) {
		let A0 = 0;
		let B0 = 1;
		let A1 = 1;
		let B1 = 1;
		let M9 = 0;
		let A2 = 0;
		let C9;
		while (Math.abs((A1 - A2)/A1) > .00001) {
			A2 = A1;
			C9 = -(A+M9)*(A+B+M9)*X/(A+2*M9)/(A+2*M9+1);
			A0 = A1+C9*A0;
			B0 = B1+C9*B0;
			M9 = M9+1;
			C9 = M9*(B-M9)*X/(A+2*M9-1)/(A+2*M9);
			A1 = A0+C9*A1;
			B1 = B0+C9*B1;
			A0 = A0/B1;
			B0 = B0/B1;
			A1 = A1/B1;
			B1 = 1;
		}
		return A1/A;
	}

	cumulativeProbability(X, df) {
		if (!df) {
			df = this.df;
		}
        let tcdf;
        if (df <= 0) {
            console.error("Degrees of freedom must be positive");
        } else {
            const A = df/2;
            const S = A + .5;
            const Z = df/(df + X*X);
            const BT = Math.exp(this.LogGamma(S) - this.LogGamma(.5) - 
                this.LogGamma(A) + A*Math.log(Z) + .5*Math.log(1 - Z));
            let betacdf;
			if (Z < (A + 1)/(S + 2)) {
                betacdf = BT*this.Betinc(Z, A, .5);
			} else {
                betacdf = 1 - BT*this.Betinc(1 - Z, .5, A);
            }
            if (X < 0) {
                tcdf = betacdf/2;
            } else {
                tcdf = 1 - betacdf/2;
            }
        }
        tcdf = Math.round(tcdf*100000)/100000;
    	return tcdf;
	}

	invCumulativeProbability(p, df) {
		if (!df){
			df = this.df;
		}
		const delta = 0.005;
        
        if (p >= 0.5) {
            let Z1 = 0;
            for (let Z = 0; Z < 100; Z++) {
                if(this.cumulativeProbability(Z, df) >= p) {
                    break;
                }
                Z1 = Z;
            }
            let Z2 = Z1;
            for (let Z = 0.0; Z < 100.0; Z += 1.0) {
                if (this.cumulativeProbability(Z1 + Z / 100.0) >= p) {
                    break;
                }
                Z2 = Z1 + (Z)/100.0;
            }
            let Z3 = Z2;
            for (let Z = 0.0; Z < 100.0; Z+=1.0) {
                if (this.cumulativeProbability(Z2 + Z / 10000.0) >= p) {
                    break;
                }
                Z3 = Z2 + (Z)/10000.0;
            }
            return Z3;
        } else {
            let Z1 = 0;
            for (let Z = 0; Z < 100; Z++) {
                if (this.cumulativeProbability(-Z, df) <= p) {
                    break;
                }
                Z1 = Z;
            }
            let Z2 = Z1;
            for (let Z = 0.0; Z < 100.0; Z += 1.0) {
                if (this.cumulativeProbability(-Z1 - Z / 100.0) <= p) {
                    break;
                }
                Z2 = Z1 + (Z) / 100.0;
            }
            let Z3 = Z2;
            for (let Z = 0.0; Z < 100.0; Z+=1.0) {
                if (this.cumulativeProbability(-Z2 - Z / 10000.0) <= p) {
                    break;
                }
                Z3 = Z2 + (Z)/10000.0;
            }
            return -Z3;
        }
    }
}

export class FDistribution {

    constructor(df1, df2) {
        this.df1 = df1;
        this.df2 = df2;
        this.EPSILON = 0.0000000001;
    }

    L504(a, f, b, iv) {
        const q = a * f / (a * f + b);
        const sa = Math.sqrt(q);
        const sl = Math.log(sa);
        const ca = Math.sqrt(1 - q);
        const cl = Math.log(ca);
        const al = Math.atan(sa / Math.sqrt(-sa * sa + 1));
        let fp = 1 - 2 * al / Math.PI;
        let r = 0.0;
        if (b != 1) {
            const c = Math.log(2 * sa / Math.PI);
            fp -= Math.exp(c + cl);
            if (b != 3) {
                const n = Math.floor((b - 3) / 2);
                for (let i = 1; i <= n; i++) {
                    const x = 2 * i + 1;
                    r += Math.log((x - 1) / x);
                    const rr = r + cl * x + c;
                    if (rr > -78.4) {
                        fp -= Math.exp(rr);
                    }
                }
            }
        }

        if (a != 1) {
            let c = r;

            if (b > 1) {
                c += Math.log(b - 1);
            }

            c += Math.log(2 / Math.PI) + sl + cl * b;

            if (c > -78.4) { fp += Math.exp(c); }

            if (a != 3) {
                const n = Math.floor((a - 3) / 2);
                r = 0;
                for (let i = 1; i <= n; i++) {
                    const x = i * 2 + 1;
                    r += Math.log((b + x - 2) / x);
                    const rr = r + sl * (x - 1) + c;
                    if (rr > -78.4) { fp += Math.exp(rr); }
                }
            }
        }
        return fp;
    }

    L401(a, f, b, iv) {
        const q = a * f / (a * f + b);
        const ql = Math.log(q);
        let fp = 0.0;
        const c = Math.log(1 - q) * b / 2;
        if (c > -78.4) {
            fp = Math.exp(c);
        }

        if (a != 2) {
            const n = Math.floor(a / 2 - 1);
            let r = 0.0;
            for (let i = 1; i <= n; i++) {
                const x = 2 * i;
                r += Math.log(b + x - 2) - Math.log(x) + ql;
                if (r + c > -78.4) {
                    fp += Math.exp(r + c);
                }
            }
        }

        if (iv == 1) {
            fp = 1 - fp;
        }

        return fp;
    }

    ProbF(dn, dd, fr) {
        let f = fr;
        let a = dn;
        let b = dd;
        let iv = 0;

        if (Math.floor(a / 2) * 2 == a) {
            //even numerator df
            const fp = this.L401(a, f, b, iv);
            return fp;
        } else if (Math.floor(b / 2) * 2 != b) {
            const fp = this.L504(a, f, b, iv);
            return fp;
        }

        f = 1 / f;
        a = dd;
        b = dn;
        iv = 1;
        return this.L401(a, f, b, iv);
    }

    cumulativeProbability(F) {
        if (this.df1 > .01 && this.df2 > .01 && F > this.EPSILON) {
            const p = 1 - this.ProbF(this.df1, this.df2, F);
            return p;
        } else{
            console.error("df1, df2, and F must be numbers greater than 0.");
        }
    }
}

export class ChiSquareDistribution {

    constructor(df) {
        this.df = df;
    }

    ChiSquaredProbability(x) {
        let a, y = 0, s, e, c, z, val;
        const df = this.df;
        const bigx = 20.0;
        const logSqrtPi = Math.log(Math.sqrt(Math.PI));
        const rezSqrtPi = 1 / Math.sqrt(Math.PI);
        if (x <= 0 || df < 1)
            return (1);
        a = 0.5 * x;
        const even = ((parseInt(2 * (df / 2), 2)) == df);
        if (df > 1)
            y = Math.exp(-a); //((-a < -bigx) ? 0.0 : Math.exp (-a));
        s = (even ? y : (2.0 * (new NormalDistribution(0.0, 1.0).cumulativeProbability(-Math.sqrt(x)))));
        if (df > 2) {
            x = 0.5 * (df - 1.0);
            z = (even ? 1.0 : 0.5);
            if (a > bigx) {
                e = (even ? 0.0 : logSqrtPi);
                c = Math.log(a);
                while (z <= x) {
                    e = Math.log(z) + e;
                    val = c * z - a - e;
                    s += Math.exp(val); //((val < -bigx) ? 0.0 : Math.exp (val));
                    z += 1.0;
                }
                return (s);
            } else {
                e = (even ? 1.0 : (rezSqrtPi / Math.sqrt(a)));
                c = 0.0;
                while (z <= x) {
                    e = e * (a / z);
                    c = c + e;
                    z += 1.0;
                }
                return (c * y + s);
            }
        } else {
            return (s);
        }
    }

    cumulativeProbability(x) {
        return 1 - this.ChiSquaredProbability(x);
    }
}
