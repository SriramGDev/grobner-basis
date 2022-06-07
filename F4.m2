Pairs = (F) -> (
	ans := {};
	for i from 0 to #F-1 do (
		for j from i+1 to #F-1 do (
			ans = append(ans, {i, j});
		);	
	);
	return ans;
);
Mon = (H) -> (
	mon := set {};
	for f in H do (
		mon = mon+set first entries monomials f;
	);
	return toList(mon);
);
SymbolicPreProcess = (G, P') -> (
	L := {};
	for p in P' do (
		f := G#(first p);
		g := G#(last p);
		L = append(L, lcm(leadMonomial(f), leadMonomial(g))//leadTerm(f)*f);
		L = append(L, lcm(leadMonomial(f), leadMonomial(g))//leadTerm(g)*g);
	);
	H := toList(set L);
	done := set apply(H, h -> leadMonomial h);
	while #(Mon(H)-done) > 0 do (
		m := max(toList(Mon(H)-done));
		done = done + set {m};
		for f in G do (
			if m%leadMonomial(f) == 0 then (
				H = toList(set H + set {(m//leadMonomial(f))*f});
			); 
		);	
	);
	sortedMonomials := reverse sort toList Mon(H);
	mList := {};
	for f in H do (
		temp := {};
		for m in sortedMonomials do (
			temp = append(temp, coefficient(m,f));	
		);
		mList = append(mList, temp);
	);
	return {matrix mList, sortedMonomials};
	
);
rows = (N, sortedMonomials) -> (
	polynomials := {};
	for i from 0 to (numgens target N)-1 do (
		p := 0;
		for j from 0 to (#sortedMonomials)-1 do (
			p = p + N_(i,j)*sortedMonomials#j;	
		);
		polynomials = append(polynomials, p);
	);
	return polynomials;
);
F4 = (F) -> (
	G := F;
	t := #G-1;
	P := Pairs(G);
	sizes := {};
	while #P != 0 do (
		P' := {};
		for p in P do (
			if random 2 == 1 then (P' = append(P', p));
		);
		if #P' == 0 then P' = append(P', first P);
		P = toList(set(P)-set(P'));
		M := SymbolicPreProcess(G, P');	
		sizes = append(sizes, {numgens target first M, numgens source first M});
		M' := reducedRowEchelonForm first M;
		R := rows(M', last M);
		LMGIdeal := ideal apply(G, f -> (if f == 0 then 0 else leadMonomial f));
		N := {};
		for f in R  do (
			if f != 0 and leadMonomial f % LMGIdeal != 0 then (
				N = append(N, f);
			);
		);
		for f in N do (
			t = t+1;
			G = append(G, f);
			for i from 0 to t-1 do (
				P = toList(set P + set {{i,t}});
			);	
		);
	);
	return {G, sizes};
);
