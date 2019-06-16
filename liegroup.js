// UNCLASSIFIED

const {Copy,Log} = require("../enum");
		
function $(N, cb) {
	var A = new Array(N);
	if (cb) for (var n=0; n<N; n++) cb( n,A );
	return A;
}

[ 
	function $(cb) {
		for (var n=0, N=this.length; n<N; n++) cb( n, this );
		return this;
	}
].Extend(Array);
	
var LG = module.exports = {
	group: GROUP,
	
	// Provide imaging methods
	
	image: function(X, M) {
	// Return supplied KxK image string X centered in an image A of M^2 entries..
	
		var A = $(M*M, (n,A) => A[n] = 0);

		for (var n=1,K=0, N=X.length; n<N; n++,K++)
			if ( X.charAt(n) == "\n" ) break;
		
		var pad = (M-K)/2;	// left, right, top, bottom padding
		
		for (var n=1, N=X.length, m=pad*M+pad; n<N; n++) {
			if ( (x = X.charAt(n)) == "\n" ) 
				m += 2*pad;
			
			else 
				A[m++] = parseInt(x);
		}

		//Log({K: K, pad: pad, image: A});
		
		return A;
	},
	
	pairs: function (M, rho) {
	/*
	Pair (xm,ym) by reflecting xm about the rho-rotational symmetry lying in an M=N^2 image.
	*/
		const {round, floor,PI,cos,sin,sqrt} = Math;
		
		var 
			N = round( sqrt(M) ),	// image A has M = NxN elements
			M2 = floor((M-N)/2)+N, 		// number of pairs in image A
			N2 = (N-1)/2,

			alpha = cos(rho * PI/180), // line of reflection symmetry
			beta = sin(rho * PI/180);
		
		return $(M, (n, P) => {
			var
				xm = n,
				x0 = [ floor(xm / N), xm % N ],  // start [image row col]
				x = [ x0[0]-N2, x0[1]-N2 ], 	// image location

				t = alpha * x[0] + beta * x[1],  // parametric parms
				u = beta * x[0] - alpha * x[1],

				y = [ alpha*t - beta*u, beta*t + alpha*u ],  // end image location
				y0 = [ round(y[0]+N2), round(y[1]+N2) ],

				ym = y0[0] * N + y0[1];
			
			P[n] = {x: xm, y: ym, x0: x0, y0: y0};
		});
	},
						
	scatter: function (A,rho,leg,cb) {
	/*
	Compute scattering of an image A about the rho-defined reflection symmetry from the named leg.  Each pair is 
	passed to the callback cb(pair) to compute its (+/-) scattering coefficients.  After all
	pairs are made, the computed scattering coefficients coff are passed to cb(coef, leg).
	*/
		const {round, floor, sqrt} = Math;
		
		var 
			M = A.length, N = round( sqrt(M) ),	// image A is NxN
			M2 = floor((M-N)/2)+N, 		// number of pairs in image A
			pairs = LG.pairs(M, rho),	// (x,y) pairs
			pos = "+", neg = "-",	// symbols +/- = sum/dif computed on this leg
			coef = { n: 0, map: $(M) };
			
		coef[pos] = $(M2);
		coef[neg] = $(M2);
		//Log( [rho, [M,N],[M2,N2],[alpha,beta] ] );

		pairs.forEach( pair => { // generate (x,y) image pairs from (x,y) pairs
			var
				img = { x: A[pair.x], y: A[pair.y] };
			
			if ( !coef.map[pair.x] ) {  // pair has not yet been connected
				coef.map[pair.x] = pair.y;  // map pair being connected
				coef.map[pair.y] = pair.x;

				if ( cb ) {  // cb computes the pair scattering symmetries (x,y)
					cb(img);

					coef[pos][coef.n] = img[pos];
					coef[neg][coef.n] = img[neg];
					coef.n++;
				}

				//Log([ x0, y0, coef.n, img ]);
			}
		});
		
		if (cb) {  // return symmetries to callback 
			cb( coef[pos] , leg + pos );
			cb( coef[neg] , leg + neg );
		}
	},
	
	haar: function haar(A,rho,depth,leg,cb) {
	/*
	Deep haar scatter image A about the rho-defined reflection symmetry to the 
	requested depth starting from the named leg.
	*/
			
		function recurse(A,depth,leg) { // pad image A to square and pass to haar 
		/* 
		At level depth, the M = NxN image A is split into two halves of H = (M-N)/2 elements
		each.  At the next depth-1 level, the KxK image is padded so that it contains no 
		more than H elements, i.e. so M+pad = KxK.
		*/
			const {round, sqrt, max, min} = Math;
			
			var
				M = A.length,
				pad = max(0, round( sqrt(M) )**2 - M),
				pads = $( pad, (n,p) => p[n] = 0 );

			//Log(['pad',depth,leg,M,pad,M+pad]);

			haar(A.concat(pads), rho, depth, leg, cb);
		}

		const {abs} = Math;
		
		LG.scatter(A, rho, leg, function (pair, leg) { 	// get scattering symmetries

			if (leg) {  // (image,leg) provided so recurse down haar tree
				//Log("depth", depth, pair.length);

				if (depth)	// recurse to next depth
					recurse( pair , depth-1, leg );

				else  	{ // callback with scatterings
					cb( pair , leg );
				}
			}
			
			else {  // (pair) provided so compute its haar scattering
				pair['+'] = pair.x + pair.y;			// sum (+)
				pair['-'] = abs( pair.x - pair.y );	// dif (-)
			}

		});
	}
}

function GROUP(N) {	// N-point group generator
/*
	Generate the 2N symmetries of an N-point group G.  For every g in G
 
		G[g] = point permutation [1:N]

	where the permutators H of G

		H[g] = op( arg[g] , point perm )

	have arguments

		arg[g] = 0,1, ... K

	There are, for N = 4 points, 8 = 2N elements in the G(4) group: 3 rotators
	(r1, r2, r3) = (45, 90, 135 degs) = rho^[1,2,3], 2 flips (f0, f1) = (H, V), 2 mirrors 
	(m0, m1) = (/, \) and an identity (e).  For N=4, arg[g] <= 3 for all g in G(4).  The 
	G(3) group contains 6 = 2N elements: 2 rotators (r1, r2) = (60, 120 degs), 3 
	swaps	(s0, s1, s2), and the identity (e).  Note odd-point groups contain swaps 
	vs the flips and mirrors of	even-point groups. P(4) confirms, for example, that a 
	rotation and flip:
	
		r1 * f0 = m0
		r1 * f1 = m1
		
	implements the mirrors, whereas two flips
	
		f0 * f1 = r2 = rho^2
		
	implements the rho^2=90 deg rotation.

	Also generated are the products P(N) for f,g,h in group G(N):

		P[ f * g ] = h

	inverses I

		f * I[ f ] = e

	involutes V

		V[ f ] * V[f] = e

	and equivalency-tests X

		X[ h ][ f*g ] = true if f*g = h	
		
	Conjugacy classes C are also computed, these involving the triple products
	_g * f * g = h, for all f,g != e in G where _g = inv(g).  We say that f and h are in the 
	same conjugacy class k when C[k][f] = C[k][h] = k.	
	*/
	 
	function rot(i,x) {  // rotation permutation
		return $(N, (n,rtn) => rtn[n] = x[ (N-i+n) % N ] );
	}

	function mirror(i,x) {  // mirror permutation
		var rtn = $(N);
		rtn[i] = x[i];
		i += N2;
		rtn[i] = x[i];
		for (var n=0, iL=i-1, iR=(i+1)%N, N1=N2-1; n<N1; n++,iL--,iR=(iR+1)%N) {
			rtn[iL] = x[iR];
			rtn[iR] = x[iL];
		}
		return rtn;
	}

	function swap(i,x) { // swap permutation
		var rtn = $(N);
		rtn[i] = x[i];
		for (var n=1,iL=(N+i-n)%N,iR=(i+n)%N; n<=N2; n++,iL=(N+i-n)%N,iR=(i+n)%N) {
			rtn[iL] = x[iR];
			rtn[iR] = x[iL];
		}
		return rtn;
	}

	function flip(i,x) { // flip permuation
		var rtn = $(N);
		for (var n=0,iL=i,iR=i+1,iN=i+1; n<iN; n++,iL--,iR++) {
			rtn[iL] = x[iR];
			rtn[iR] = x[iL];
		}
		i += N2;
		for (var n=0,iL=i,iR=i+1,iN=N-i-1; n<iN; n++,iL--,iR++) {
			rtn[iL] = x[iR];
			rtn[iR] = x[iL];
		}
		return rtn;
	}

	function ident(i,x) { // identity permutation
		return x;
	}

	function find(x, cb) { // callback cb(h) with element h of G that produces perm x
		function eq(x,y) { // test if permuations x,y are equal
			for (var n=0;n<N;n++)
				if (x[n] != y[n]) return false;

			return true;
		}

		for (var h in G)
			if ( eq(x, G[h]) )
				return cb(h);

		Log("Houston we have a problem - G is not a group!");
	}

	function index(k) {
		k = k || 1;
		return $(N, (n,rtn) => rtn[n] = n*k );
	}

	var 
		G = this.G = {},
		H = this.H = {},
		P = this.P = {},
		X = this.X = {},
		I = this.I = {},
		A = this.A = {},
		V = this.V = [],
		C = this.C = [{e:0}],
		odd = this.odd = (N % 2) ? true : false,
		even = this.even = ! this.odd,
		e = G.e = [],
		rho =  this.rho = index(180/N),
		order = this.order = 2*N,
		sym = this.sym = { // symmetry labels
			std: {},
			fav: {
				f0: "H",	// flip
				f1: "V",
				m0: "/",	// mirror
				m1: "\\", 
				r0: "N"	// rotation
			}
		},
		N2 = even ? N/2 : (N-1)/2;

	this.moves = {flips:even?N:0, mirrors:even?N-2:0, swaps:odd?N:0};

	for (var n=1;n<=N;n++) e.push(n); H[g="e"] = ident; A[g] = 0;  // H[g] = op(arg,perm)

	for (var n=1;n<N;n++) G[g="r"+n] = (H[g]=rot)(A[g]=n,e);  // G[g] = perm

	if (even) {
		for (var n=0; n<N2; n++) G[g="f"+n] = (H[g]=flip)(A[g]=n,e);
		for (var n=0; n<N2; n++) G[g="m"+n] = (H[g]=mirror)(A[g]=n,e);
	}

	else
		for (var n=0; n<N; n++) G[g="s"+n] = (H[g]=swap)(A[g]=n,e);

	for (var f in G) for (var g in G)
		find( fg = H[g](A[g], G[f]), h => {  // H[g] = PermOp(arg, perm)
			if (h == "e") {					// if h = identity returned
				I[f] = g; I[g] = f;				// save inverses  f * I[ f ] = e
				if (f == g) V.push(f);		// save involutes V[ f ] * V[ f ] = e
			}

			P[fg = f+"*"+g] = h; 	// save products f*g = h

			if (! X[h] ) X[h] = {}; 	// reserve for tests 
			X[h][fg] = (f[0] != "r" && g[0] != "r") ? true : false;  // save X[ h ][ f*g ] true if f*g = h
		});  

	// Generate conjugacy classes C.

	for (var f in G) for (var g in G) if (f != "e" && g != "e") {
		var _g = I[ g ], _gf = P[ _g+"*"+f ], _gfg = H[g](A[g], G[_gf]);

		find( _gfg, h => {
			for (var k=0, K=C.length; k<K; k++)
				if ( C[k][f] ) return C[k][h] = k;
				//else
				//if (C[k][h] ) return C[k][f]=k;

			C.push( {} );
			C[K][f] = C[K][h] = K+1;
		});
	}
}

switch ( process.argv[2] ) { //< unit tests
	case "L3":
		Log("G(4) has 8 elements", new GROUP(4));
		//Log("G(3) has 6 elements", new GROUP(3));
		//Log("G(5) has 10 elements", new GROUP(5));
		break;
		
	case "L2":
		//Log( "4x4 V=f0 sym pairs", LG.pairs(16, 0 ) );
		//Log( "4x4 H=f1 sym pairs", LG.pairs(16, 90 ) );
		//Log( "4x4 /=m0 sym pairs", LG.pairs(16, 45 ) );
		Log( "4x4 \\=m1 sym pairs", LG.pairs(16, -45 ) );
		break;
		
	case "L1":
	/*7811
	6921
	5431
	2222
	*/
		var 
			N = 4,
			depth = 3,
			G = new GROUP(N),
			Uset = {},
			A = LG.image(`
1789
1234
4321
6543
`, 16);		// 256=16x16 image from string

		//Log(A);
		//Log(G);
		
		LG.haar( A,	G.rho[1], depth , "", function (S,leg) {

			function dot(a,b) {
				var
					sum = 0;
				
				for (var n=0, N=a.length; n<N; n++) sum += a[n] * b[n];
				return sum;
			}

			function proj(u,v) {
				return scale( copy(u), - dot(v, u) / dot(u, u) );
			}

			function scale(u,a) {
				return u.$( (n,u) => u[n] *= a );
			}

			function add(u,v) {
				return u.$( (n,u) => u[n] += v[n] );
			}

			function copy(u) {
				return $(u.length, (n,x) => x[n] = u[n] );
			}

			function gs(v, uset) {
				function allZero(u) {	
					const {abs} = Math;

					for (var n=0, N=u.length; n<N; n++)
						if ( abs(u[n]) > 1e-3 ) 
							return false;

					return true;
				}

				var  u = copy(v);

				for (var leg in uset) 
					if ( allZero( add(u, proj(uset[leg], v)) ) ) {
						//Log(["drop "+n, v]);
						return null;
					}

				return u;
			}

			Log(leg, S.length);

			if (false) 
				Uset[leg] = copy(S);

			else
			if (u = gs(S, Uset)) 
				Uset[leg] = copy(u);

		});

		//Log(Uset);
		for (var n in Uset) Log(n);
/* 
Here, this produces 5 significant scattering vectors of length 28 in 
each leg of length 4 = depth+1:

		++++
		+++-
		++-+
		+-++
		-+++

where +/- denotes sum/dif scattering calculations along each leg of the scattering.  
Here the M = 256 = 16^2 image at depth 4 partitions into halves of length 136, 78, 45, 28 
at depth = 3, 2, 1, 0; and thus images of length 144 = 12^2 = 136 + 8, 81 = 9^2 = 78 + 3, 
49 = 7^2 = 45 + 4, 36 = 6^2 = 28 + 8.  Thus, here, the last (depth=0) image (aka scattering
vector S) is of length 28.  Of the 16 = 2^(depth+1) legs taken, only 5 produced significant 
scattering, as measured by projecting S on the other 16 vectors.
*/
		break;
}

// UNCLASSIFIED
				   
