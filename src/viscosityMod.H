cout<<"Bad Mood\n";
if(TimeCounter == 200)
{
forAll(mesh.C(), celli)	{
	viscoDPD[celli] = shearRate[celli] * 0.01;
}
}
