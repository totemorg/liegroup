// UNCLASSIFIED

function GEN(N) {
	function rot(i) {
		var rtn = new Array(N);
		for (var n=0; n<N; n++)
			rtn[n] = e[ (N-i+n) % N ];
		return rtn;
	}
	
	function mirror(i) {
		var rtn = new Array(N);
		rtn[i] = e[i];
		i += N2;
		rtn[i] = e[i];
		for (var n=0, iL=i-1, iR=(i+1)%N, N1=N2-1; n<N1; n++,iL--,iR=(iR+1)%N) {
			rtn[iL] = e[iR];
			rtn[iR] = e[iL];
		}
		return rtn;
	}
		
	function swap(i) {
		var rtn = new Array(N);
		rtn[i] = e[i];
		for (var n=1,iL=(N+i-n)%N,iR=(i+n)%N; n<=N2; n++,iL=(N+i-n)%N,iR=(i+n)%N) {
			rtn[iL] = e[iR];
			rtn[iR] = e[iL];
		}
		return rtn;
	}

	function flip(i) {
		var rtn = new Array(N);
		for (var n=0,iL=i,iR=i+1,iN=i+1; n<iN; n++,iL--,iR++) {
			rtn[iL] = e[iR];
			rtn[iR] = e[iL];
		}
		i += N2;
		for (var n=0,iL=i,iR=i+1,iN=N-i-1; n<iN; n++,iL--,iR++) {
			rtn[iL] = e[iR];
			rtn[iR] = e[iL];
		}
		return rtn;
	}

	var 
		X = this.X = {},
		G = this.G = {},
		odd = this.odd = (N % 2) ? true : false,
		even = this.even = ! this.odd,
		e = G.e = [],
		rho =  this.rho = 360/N,
		rho2 = rho/2,
		N2 = even ? N/2 : (N-1)/2;
	
	this.moves = {flips:even?N:0, mirros:even?N-2:0, swaps:odd?N:0};
	
	for (var n=1;n<=N;n++) e.push(n);
	for (var n=1;n<N;n++) G["r"+n] = rot(n);
	
	if (even) {
		for (var n=0; n<N2; n++) G["f"+n] = flip(n);
		for (var n=0; n<N2; n++) G["m"+n] = mirror(n);
	}
	else
		for (var n=0; n<N; n++) G["s"+n] = swap(n);
	
}

console.log(N=5, new GEN(N));

// CLASSIFIED
				   