
	while (1) {
		if (me == 0) {
			if (fgets(line,1024,sp) == NULL) n = 0;
			else n = strlen(line) + 1;
			if (n == 0) fclose(sp);
	}
		MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
		if (n == 0) break;
		MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
		if (lammps == 1) lmp->input->one(line);
	}
