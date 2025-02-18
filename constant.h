constexpr double kMe = 0.510998950; // MeV
constexpr double kMp = 938.2720813; // MeV
constexpr double kMn = 939.56542052; // MeV
const double kamu=931.49410242;//MeV/c2
const double bez=15.73e-6;//Binding energy of electrons, MeV

auto *geom = new TGeoManager("geom", "geometry");
auto *gElementTable = geom->GetElementTable();

bool getNucleusName1993(short a,short z,string &name) {
	TGeoElementRN *gElementRN = gElementTable->GetElementRN(a,z);
	if(gElementRN ==nullptr) return false;
	TString label(gElementRN->GetName());
	TObjArray* obj = label.Tokenize("-");
	name=std::string(((TObjString*)obj->At(1))->String());
	return true;
}
//calc nuclear mass, see https://web-docs.gsi.de/~wolle/TELEKOLLEG/KERN/LECTURE/Fraser/L1.pdf
bool getNucleusMass1993(short a,short z,double &m){
	TGeoElementRN *gElementRN = gElementTable->GetElementRN(a,z);
	/*if(a==64 && z==33){
		double mex=-39.530;//MeV
		m= a * kamu +mex -z*kMe +bez *pow(z, 7./3);
		return true;
	}*/
	if(gElementRN ==nullptr){
		//cout<<RED<<"No info for:\t"<<GREEN<<"(A="<<a<<", Z="<<z<<RESET<<")\n";
		return false;
	}
	//if((gElementRN->Stable())) return false;//stable?
	//if(gElementRN->HalfLife() < 1e-100) return false;
	double mex = gElementRN->MassEx();//mass excess
	m= a * kamu +mex -z*kMe +bez *pow(z, 7./3);
	return true;
}
double getNucleusMass1993(short a,short z){
	TGeoElementRN *gElementRN = gElementTable->GetElementRN(a,z);
	if(gElementRN ==nullptr){
		//cout<<RED<<"No info for:\t"<<GREEN<<"(A="<<a<<", Z="<<z<<RESET<<")\n";
		return -1.;
	}
	double mex = gElementRN->MassEx();//mass excess
	return  a * kamu +mex -z*kMe +bez *pow(z, 7./3);
}
bool getBindingEnergy1993(short a,short z,double &be){
	double m=0.;
	if(getNucleusMass1993(a,z,m)){
		be=z*kMp + (a-z)*kMn -m;
		return true;
	}
	else return false;
}
bool nucleiType(const short &z,const short &n,short &type){
	if(z%2==0 && n%2==0)//even-Z even-N nuclides
		type=0;//MeV
	else if(z%2==0 && n%2==1)//even-Z odd-N nuclides
		type=1;//MeV
	else if(z%2==1 && n%2==0)//odd-Z even-N nuclides
		type=2;//MeV
	else if(z%2==1 && n%2==1)//odd-Z odd-N nuclides
		type=3;//MeV
	else return false;
	return true;
}
