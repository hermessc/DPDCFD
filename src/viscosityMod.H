cout<<"Bad Mood\n";
if(TimeCounter == 20)
{
forAll(mesh.C(), celli)	{
	viscoDPD[celli] = shearRate[celli] * 0.01;
	storageMatrix[0][0] = 1000;
}
}
