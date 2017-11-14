forAll(mesh.C(),celli) {
	for (int i = 0; i < 100; i++) {
		if (std::round(shearRate[celli]) == bigMM[i][0]) {
		viscoDPD[celli] = bigMM[i][1];
		}

	}
}
