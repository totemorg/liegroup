// UNCLASSIFIED

function GROUP(N) {
	function rot(i,x) {  // rotation permutation
		var rtn = new Array(N);
		for (var n=0; n<N; n++)
			rtn[n] = x[ (N-i+n) % N ];
		return rtn;
	}
	
	function mirror(i,x) {  // mirror permutation
		var rtn = new Array(N);
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
		var rtn = new Array(N);
		rtn[i] = x[i];
		for (var n=1,iL=(N+i-n)%N,iR=(i+n)%N; n<=N2; n++,iL=(N+i-n)%N,iR=(i+n)%N) {
			rtn[iL] = x[iR];
			rtn[iR] = x[iL];
		}
		return rtn;
	}

	function flip(i,x) { // slip permuation
		var rtn = new Array(N);
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
	
	function eq(x,y) { // test permuations are equal
		for (var n=0;n<N;n++)
			if (x[n] != y[n]) return false;
		
		return true;
	}
	
	function find(x, cb) { // find permuation and pass to callback
		for (var h in G)
				if ( eq(x, G[h]) )
					return cb(h);
		
		console.log("Houston we have a problem - G is not a group!");
	}
	
	function index(k) {
		k = k || 1;
		var rtn = new Array(N);
		for (var n=0; n<N; n++) rtn[n] = n*k;
		return rtn;
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
				f0: "h",
				f1: "v",
				m0: "m",
				n0: "n"
			}
		},
		N2 = even ? N/2 : (N-1)/2;
	
	this.moves = {flips:even?N:0, mirrors:even?N-2:0, swaps:odd?N:0};
	
	// generate the 2N symmetries of G and their corresponding permutators H with arguments A.
	
	for (var n=1;n<=N;n++) e.push(n); H[g="e"] = ident; A[g] = 0;
	
	for (var n=1;n<N;n++) G[g="r"+n] = (H[g]=rot)(A[g]=n,e);
	
	if (even) {
		for (var n=0; n<N2; n++) G[g="f"+n] = (H[g]=flip)(A[g]=n,e);
		for (var n=0; n<N2; n++) G[g="m"+n] = (H[g]=mirror)(A[g]=n,e);
	}
	else
		for (var n=0; n<N; n++) G[g="s"+n] = (H[g]=swap)(A[g]=n,e);
	
	// generate products P, inverses I, involutes V and rotation-tests X
	
	for (var f in G) for (var g in G)
		find( fg = H[g](A[g], G[f]), function (h) {
			if (h == "e") {
				I[f] = g; I[g] = f;
				if (f == g) V.push(f);
			}
			P[fg = f+"*"+g] = h;
			if (! X[h] ) X[h] = {}; X[h][fg] = (f[0] != "r" && g[0] != "r") ? 1 : 0;
		});
	
	// generate conjugacy classes C
	
	for (var f in G) for (var g in G) if (f != "e" && g != "e") {
		var _g = I[g], _gf = P[_g+"*"+f], _gfg = H[g](A[g], G[_gf]);
		
		find(_gfg, function (h) {
			for (var k=0, K=C.length; k<K; k++)
				if ( C[k][f] ) return C[k][h]=k;
				//else
				//if (C[k][h] ) return C[k][f]=k;
			
			C.push( {} );
			C[K][f] = C[K][h] = K+1;
		});
	}
	
	var image = this.image = function(X, M) { // return KxK image X centered in NxN image
		var A = new Array(M*M);

		for (var n=0,N=A.length; n<N; n++) A[n] = 0;

		for (var n=1,K=0, N=X.length; n<N; n++,K++)
			if ( X.charAt(n) == "\n" ) break;
		
		var pad = (M-K)/2;
		
		for (var n=1, N=X.length, m=pad*M+pad; n<N; n++)
			if ( (x = X.charAt(n)) == "\n" ) 
				m += 2*pad;
			else
				A[m++] = parseInt(x);

		//console.log(A);
		return A;
	}			
		
	// Pass the computed scattering symmetries of an image A about angle rho to a callback cb.
	
	var scatter = this.scatter = function (A,rho,leg,cb) {
		var 
			M = A.length, N = Math.round( Math.sqrt(M) ),	// image A is NxN
			M2 = Math.floor((M-N)/2)+N, N2 = (N-1)/2, 		// number of pairs in image A

			c = Math.PI / 180, 	// compute line of reflection symmetry
			alpha = Math.cos(rho * c), 
			beta = Math.sin(rho * c),

			syms = { x: new Array(M2), y: new Array(M2), n: 0, map: new Array(M) };

		//console.log( [rho, [M,N],[M2,N2],[alpha,beta] ] );

		/*
		var
			a = N2 * Math.sqrt( 2 ),
			a2 = a*a,		
			dx = dy = 1,
			dt = Math.max(alpha, beta),
			gamma = alpha*alpha + beta*beta,			
			du = Math.abs( dy*alpha/gamma - dx*beta/gamma );
			
		for (var t=-a; t<a; t+=dt)
			for (var b=Math.sqrt(a2 - t*t), u=0; u<b; u+=du) {
				
				var 
					x = [ alpha*t - beta*u, beta*t + alpha*u ], // coordinates of pair
					y = [ alpha*t + beta*u, beta*t - alpha*u ],
					
					x0 = [ Math.floor(x[0] + N2), Math.floor(x[1] + N2) ],  // shift to center of image
					y0 = [ Math.floor(y[0] + N2), Math.floor(y[1] + N2) ],
					
					xm = N * x0[0] + x0[1],  // index into image
					ym = N * y0[0] + y0[1],
					
					pair = { x: A[xm], y: A[ym] };

				console.log( [ x,x0,xm, y,y0,ym, pair ]);

				if ( pair.x != undefined && pair.y != undefined ) { // pair inside image
					if (cb) cb(pair);

					A[y0][y1] = pair.x;
					A[x0][x1] = pair.y;
					//console.log([x0,x1], [y0,y1], pair);
					//console.log([ x,y, pair]);
				}
				//else
				//	console.log([x0,x1], [y0,y1]);
					
			}
		*/

		for (var xm=0; xm<M; xm++) {  // generate (x,y) pairs
			var
				x0 = [ Math.floor(xm / N), xm % N ],  // pair start [image row col]
				x = [ x0[0]-N2, x0[1]-N2 ], 	// image location

				t = alpha * x[0] + beta * x[1],  // pair parametric parms
				u = beta * x[0] - alpha * x[1],

				y = [ alpha*t - beta*u, beta*t + alpha*u ],  // pair end image location
				y0 = [ Math.round(y[0]+N2), Math.round(y[1]+N2) ],
				ym = y0[0] * N + y0[1],

				pair = { x: A[xm], y: A[ym] };
			
			if ( !syms.map[xm] ) {  // pair has not been connected
				syms.map[xm] = ym;  // map pair being connected
				syms.map[ym] = xm;

				if ( cb ) {  // cb computes the pair scattering symmetries (x,y)
					cb(pair);

					syms.x[syms.n] = pair.x;
					syms.y[syms.n] = pair.y;
					syms.n++;
				}

				//console.log([ x0, y0, syms.n, [pair.x, pair.y] ]);
			}
		}
		
		if (cb) {  // return symmetries to callback 
			cb( syms.x , "x");
			cb( syms.y , "y");
		}
	}
	
	// deep haar scatter image A about rho symmetry until depth level
	
	var haar = this.haar = function (A,rho,level,leg,cb) {
		
		function recurse(A,level,leg) { // pad image A to square and pass to haar
			var
				M = A.length,
				pad = Math.max(0, Math.pow(Math.round( Math.sqrt(M) ),2) - M),
				pads = new Array(pad);

			//console.log(['pad',level,leg,M,pad]);
			for (var n=0; n<pad; n++) pads[n] = 0;

			haar(A.concat(pads), rho, level, leg, cb);
		}

		scatter(A, rho, leg, function (pair, lab) { 	// get scattering symmetries

			if (pair.constructor == Array)   // image so recurse down haar tree
				if (level)  					// recurse to next level
					recurse( pair , level-1, leg+lab);

				else  							// callback with scatterings
					cb( pair , leg+lab );
			
			else {  // pair so compute its haar scattering
				var x = pair.x, y = pair.y;
				
				pair.x = x+y;
				pair.y = Math.abs(x-y);
			}

		});
	}
	
}

var G = new GROUP(N=4);

/*7811
6921
5431
2222
*/
var A = G.image(`
1111
1111
1111
1111
`, 16);

var NS = 0;

G.haar( A,	G.rho[1], 3 , "", function (S,leg) {
	console.log([NS++, leg, S]);
});

// UNCLASSIFIED
				   