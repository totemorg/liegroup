// UNCLASSIFIED

function GEN(N) {
	function rot(i,x) {
		var rtn = new Array(N);
		for (var n=0; n<N; n++)
			rtn[n] = x[ (N-i+n) % N ];
		return rtn;
	}
	
	function mirror(i,x) {
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
		
	function swap(i,x) {
		var rtn = new Array(N);
		rtn[i] = x[i];
		for (var n=1,iL=(N+i-n)%N,iR=(i+n)%N; n<=N2; n++,iL=(N+i-n)%N,iR=(i+n)%N) {
			rtn[iL] = x[iR];
			rtn[iR] = x[iL];
		}
		return rtn;
	}

	function flip(i,x) {
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

	function ident(i,x) {
		return x;
	}
	
	function eq(x,y) {
		for (var n=0;n<N;n++)
			if (x[n] != y[n]) return false;
		
		return true;
	}
	
	function find(x, cb) {
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
		sym = this.sym = {
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
	
	// gen the 2N symmetries of G and their corresponding permutators H with arguments A.
	
	for (var n=1;n<=N;n++) e.push(n); H[g="e"] = ident; A[g] = 0;
	
	for (var n=1;n<N;n++) G[g="r"+n] = (H[g]=rot)(A[g]=n,e);
	
	if (even) {
		for (var n=0; n<N2; n++) G[g="f"+n] = (H[g]=flip)(A[g]=n,e);
		for (var n=0; n<N2; n++) G[g="m"+n] = (H[g]=mirror)(A[g]=n,e);
	}
	else
		for (var n=0; n<N; n++) G[g="s"+n] = (H[g]=swap)(A[g]=n,e);
	
	// gen product P, inverse I, involute V and the rotation-tested X maps
	
	for (var f in G) for (var g in G)
		find( fg = H[g](A[g], G[f]), function (h) {
			if (h == "e") {
				I[f] = g; I[g] = f;
				if (f == g) V.push(f);
			}
			P[fg = f+"*"+g] = h;
			if (! X[h] ) X[h] = {}; X[h][fg] = (f[0] != "r" && g[0] != "r") ? 1 : 0;
		});
	
	// gen conjugacy classes C
	
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
	
	// define the the reflection of an image A about the n'th mirror
	
	this.matrix = function(X) {
		var A = [], B = [];

		for (var n=1, N=X.length; n<N; n++)
			if ((x = X.charAt(n)) == "\n") {
				A.push(B);
				B = new Array();
			}
			else
				B.push( parseInt(x) );

		//console.log(A);
		return A;
	}			
		
	this.reflector = function (A,n,s) {		
		var 
			M = A.length, N = A[0].length, 
			M2 = M/2, N2 = N/2,
			a = Math.sqrt( M2*M2 + N2*N2 ) *(s||1),
			a2 = a*a,
			c = Math.PI / 180,
			line = [ alpha = Math.cos(rho[n] * c), beta = Math.sin(rho[n] * c) ],
			dt = Math.max(alpha, beta),
			dx = dy = 1,
			gamma = alpha*alpha + beta*beta,
			du = Math.abs( dy*alpha/gamma - dx*beta/gamma );
		
		console.log([[M,N], a, line, [dt,du], [alpha, beta] ]);
		
		for (var t=-a; t<a; t+=dt)
			for (var b=Math.sqrt(a2 - t*t), u=0; u<b; u+=du) {
				
				var 
					x = [ alpha*t - beta*u, beta*t + alpha*u ],
					y = [ alpha*t + beta*u, beta*t - alpha*u ],
					
					x0 = Math.trunc(x[0] + M2), x1 = Math.trunc(x[1] + M2),
					y0 = Math.trunc(y[0] + N2), y1 = Math.trunc(y[1] + N2),
				
					yHold = (A[y0] || [])[y1],
					xHold = (A[x0] || [])[x1];
								
				if ( xHold != undefined && yHold != undefined ) {
					A[y0][y1] = xHold;
					A[x0][x1] = yHold;
					//console.log([x0,x1],[y0,y1],[xHold,yHold]);
				}
				else
					console.log([x0,x1],[y0,y1]);
					
			}
		
		return A;
	}
	
	this.test = `
00000
07810
06920
05430
00000
`;
		
}

var p = new GEN(N=4);

//console.log(p);
var A = p.matrix(p.test);
console.log(p.reflector( A, 1) );

// CLASSIFIED
				   