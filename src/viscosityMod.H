//H
// modificatore di viscosità a scopo dimostrativo fittizio, va rimpiazzato con una legge che converta le unità dpd in unità fisiche

cout<<"Bad Mood\n";
if(TimeCounter == 200)
{
forAll(mesh.C(), celli)	{
	viscoDPD[celli] = shearRate[celli] * 0.01;
	//storageMatrix[0][0] = 1000;
}
}
