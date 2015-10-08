enum genealsourcecode
{
	MUON, TPA, TPB, CALOMUON, TPC // general TPs
};

enum sourcecode
{
	CALO, MUONS, MUID, // muon...
	MUONSTPA, MUONSTPB, MUONSTPC, // muons chain TPs
	MUIDTPA, MUIDTPB, MUIDTPC // muid chain TPs
};

enum categorycodes
{
	MUONS3,
	MUONS2TPA1, MUONS2TPB1, MUONS2CALO1,
	MUONS1TPA2, MUONS1TPB2,
	MUONS1TPA1TPB1,
	MUONS0TPA3, MUONS0TPB3,
	MUONS0TPA2TPB1, MUONS0TPA2CALO1,
	MUONS0TPA1TPB2,
	
	MUID3,
	MUID2TPA1, MUID2TPB1, MUID2CALO1,
	MUID1TPA2, MUID1TPB2,
	MUID1TPA1TPB1,
	MUID0TPA3, MUID0TPB3,
	MUID0TPA2TPB1, MUID0TPA2CALO1,
	MUID0TPA1TPB2
};


void buildTPmus()
{	
	tpmus.insert(make_pair(CALO,     "CaloMuonCollection"));
	tpmus.insert(make_pair(MUID,     "MuidMuonCollection"));
	tpmus.insert(make_pair(MUONS,    "Muons"));
	tpmus.insert(make_pair(MUONSTPA, "CombinedFitMuonParticles"));
	tpmus.insert(make_pair(MUONSTPB, "SegmentTagTrackParticles"));
	tpmus.insert(make_pair(MUONSTPC, "MuGirlRefittedTrackParticles"));
	tpmus.insert(make_pair(MUIDTPA,  "MuidCombTrackParticles"));
	tpmus.insert(make_pair(MUIDTPB,  "MuTagIMOTrackParticles"));
	tpmus.insert(make_pair(MUIDTPC,  "MuGirlRefittedTrackParticles")); // same as for MUONS
}


TString classifyTripletShort(unsigned int vtx)
{
	TString type = classifyTriplet(vtx);
	TString shorttype = "";
	if     (type==tpmus[MUONS]+":"+tpmus[MUONS]+":"+tpmus[MUONS])          shorttype = "3muons";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONS]+":"+tpmus[MUONSTPA])       shorttype = "2muons1tpmuA";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONS]+":"+tpmus[MUONSTPB])       shorttype = "2muons1tpmuB";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPA])    shorttype = "1muons2tpmuA";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPB])    shorttype = "1muons2tpmuB";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPB])    shorttype = "1muons1tpmuA1tpmuB";
	else if(type==tpmus[MUONS]+":"+tpmus[MUONS]+":"+tpmus[CALO])           shorttype = "2muons1calomu";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]) shorttype = "0muons3tpmuA";
	else if(type==tpmus[MUONSTPB]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPB]) shorttype = "0muons3tpmuB";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]+":"+tpmus[MUONSTPB]) shorttype = "0muons2tpmuA1tmpmuB";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPB]+":"+tpmus[MUONSTPB]) shorttype = "0muons1tpmuA2tmpmuB";
	else if(type==tpmus[MUONSTPA]+":"+tpmus[MUONSTPA]+":"+tpmus[CALO])     shorttype = "0muons2tpmuA1calomu";
	
	else if(type==tpmus[MUID]+":"+tpmus[MUID]+":"+tpmus[MUID])             shorttype = "3muid";
	else if(type==tpmus[MUID]+":"+tpmus[MUID]+":"+tpmus[MUIDTPA])          shorttype = "2muid1tpmuA";
	else if(type==tpmus[MUID]+":"+tpmus[MUID]+":"+tpmus[MUIDTPB])          shorttype = "2muid1tpmuB";
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPA])       shorttype = "1muid2tpmuA";
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPB])       shorttype = "1muid2tpmuB";
	else if(type==tpmus[MUID]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPB])       shorttype = "1muid1tpmuA1tpmuB";
	else if(type==tpmus[MUID]+":"+tpmus[MUID]+":"+tpmus[CALO])             shorttype = "2muid1calomu";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPA])    shorttype = "0muid3tpmuA";
	else if(type==tpmus[MUIDTPB]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPB])    shorttype = "0muid3tpmuB";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]+":"+tpmus[MUIDTPB])    shorttype = "0muid2tpmuA1tmpmuB";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPB]+":"+tpmus[MUIDTPB])    shorttype = "0muid1tpmuA2tmpmuB";
	else if(type==tpmus[MUIDTPA]+":"+tpmus[MUIDTPA]+":"+tpmus[CALO])       shorttype = "0muid2tpmuA1calomu";

	else _FATAL("Triplet cannot be classified");

	return shorttype;
}
int classifyTripletCode(TString shorttype)
{
	if     (shorttype=="3muons")              return MUONS3;
	else if(shorttype=="2muons1tpmuA")        return MUONS2TPA1;
	else if(shorttype=="2muons1tpmuB")        return MUONS2TPB1;
	else if(shorttype=="2muons1calomu")       return MUONS2CALO1;
	else if(shorttype=="1muons2tpmuA")        return MUONS1TPA2;
	else if(shorttype=="1muons2tpmuB")        return MUONS1TPB2;
	else if(shorttype=="1muons1tpmuA1tpmuB")  return MUONS1TPA1TPB1;
	else if(shorttype=="0muons3tpmuA")        return MUONS0TPA3;
	else if(shorttype=="0muons3tpmuB")        return MUONS0TPB3;
	else if(shorttype=="0muons2tpmuA1tmpmuB") return MUONS0TPA2TPB1;
	else if(shorttype=="0muons2tpmuA1calomu") return MUONS0TPA2CALO1;
	else if(shorttype=="0muons1tpmuA2tmpmuB") return MUONS0TPA1TPB2;
	
	else if(shorttype=="3muid")              return MUID3;
	else if(shorttype=="2muid1tpmuA")        return MUID2TPA1;
	else if(shorttype=="2muid1tpmuB")        return MUID2TPB1;
	else if(shorttype=="2muid1calomu")       return MUID2CALO1;
	else if(shorttype=="1muid2tpmuA")        return MUID1TPA2;
	else if(shorttype=="1muid2tpmuB")        return MUID1TPB2;
	else if(shorttype=="1muid1tpmuA1tpmuB")  return MUID1TPA1TPB1;
	else if(shorttype=="0muid3tpmuA")        return MUID0TPA3;
	else if(shorttype=="0muid3tpmuB")        return MUID0TPB3;
	else if(shorttype=="0muid2tpmuA1tmpmuB") return MUID0TPA2TPB1;
	else if(shorttype=="0muid2tpmuA1calomu") return MUID0TPA2CALO1;
	else if(shorttype=="0muid1tpmuA2tmpmuB") return MUID0TPA1TPB2;
	
	else _FATAL("Triplet cannot be classified");
	
	return -999;
}
int getSrcType(unsigned int vtx, unsigned int index)
{
	TString srcName = vtx_srcName->at(vtx)[index]; // this is not necessarily ordered
	int srcType = -999;
	for(TMapiTS::iterator it=tpmus.begin() ; it!=tpmus.end() ; ++it)
	{
		if(it->second==srcName) { srcType = it->first; break; }
	}
	if(srcType==-999) _FATAL("Source type could not be found");
	return srcType;
}
bool isMuonSrc(int srcType)
{
	return (srcType<=MUID);
}
bool isCaloSrc(int srcType)
{
	return (srcType==CALO);
}
bool isTPmuSrc(int srcType)
{
	return (srcType>MUID);
}
bool isTPaSrc(int srcType)
{
	return (srcType==MUONSTPA || srcType==MUIDTPA);
}
bool isTPbSrc(int srcType)
{
	return (srcType==MUONSTPB || srcType==MUIDTPB);
}
int getSrcCode(int srcType)
{
	if(isMuonSrc(srcType)) return MUON;
	if(isCaloSrc(srcType)) return CALOMUON;
	if(isTPaSrc(srcType))  return TPA;
	if(isTPbSrc(srcType))  return TPB;
	return -1;
}

void getSrc(unsigned int vtx, sources& src)
{	
	src.vtx       = vtx;
	src.type      = classifyTriplet(vtx);
	src.shorttype = classifyTripletShort(vtx);
	for(int i=0 ; i<3 ; ++i) src.srcName[i]  = vtx_srcName->at(vtx)[i];
	for(int i=0 ; i<3 ; ++i) src.trkIndex[i] = vtx_trkIndex->at(vtx)[i];
	for(int i=0 ; i<3 ; ++i) src.srcIndex[i] = vtx_srcIndex->at(vtx)[i];
	for(int i=0 ; i<3 ; ++i) src.srcType[i]  = getSrcType(vtx,i);
	for(int i=0 ; i<3 ; ++i) src.isMuon[i]   = isMuonSrc(src.srcType[i]);
	for(int i=0 ; i<3 ; ++i) src.isCalo[i]   = isCaloSrc(src.srcType[i]);
	for(int i=0 ; i<3 ; ++i) src.isTPmu[i]   = isTPmuSrc(src.srcType[i]);
	for(int i=0 ; i<3 ; ++i) src.isTPa[i]    = isTPaSrc(src.srcType[i]);
	for(int i=0 ; i<3 ; ++i) src.isTPb[i]    = isTPbSrc(src.srcType[i]);
	for(int i=0 ; i<3 ; ++i) src.srcCode[i]  = getSrcCode(src.srcType[i]);
	
	double px1 = vtx_reftrks_px->at(vtx)[0];
	double px2 = vtx_reftrks_px->at(vtx)[1];
	double px3 = vtx_reftrks_px->at(vtx)[2];
	double py1 = vtx_reftrks_py->at(vtx)[0];
	double py2 = vtx_reftrks_py->at(vtx)[1];
	double py3 = vtx_reftrks_py->at(vtx)[2];
	double pz1 = vtx_reftrks_pz->at(vtx)[0];
	double pz2 = vtx_reftrks_pz->at(vtx)[1];
	double pz3 = vtx_reftrks_pz->at(vtx)[2];
	TLorentzVector p1,p2,p3;
	TMapdi srcMap;
	p1.SetXYZM(px1,py1,pz1,muonMassMeV); src.srcTlv[0] = p1; double pT1 = p1.Pt(); srcMap.insert(make_pair(pT1,0));
	p2.SetXYZM(px2,py2,pz2,muonMassMeV); src.srcTlv[1] = p2; double pT2 = p2.Pt(); srcMap.insert(make_pair(pT2,1));
	p3.SetXYZM(px3,py3,pz3,muonMassMeV); src.srcTlv[2] = p3; double pT3 = p3.Pt(); srcMap.insert(make_pair(pT3,2));
	TMapdi::reverse_iterator rit=srcMap.rbegin();
	src.srcOrder[rit->second] = 1; ++rit;
	src.srcOrder[rit->second] = 2; ++rit;
	src.srcOrder[rit->second] = 3;
	
	src.p3body  = p1+p2+p3;
	src.p2body12 = p1+p2;
	src.p2body13 = p1+p3;
	src.p2body23 = p2+p3;
	
	// double q1 = (src.isMuon[0] || src.isCalo[0]) ? qtrk(trks_qoverp->at(src.trkIndex[0])) : qtrk(tpmu_vd[src1+"_qOverP"]->at(src.srcIndex[0]));
	// double q2 = (src.isMuon[1] || src.isCalo[1]) ? qtrk(trks_qoverp->at(src.trkIndex[1])) : qtrk(tpmu_vd[src2+"_qOverP"]->at(src.srcIndex[1]));
	// double q3 = (src.isMuon[2] || src.isCalo[2]) ? qtrk(trks_qoverp->at(src.trkIndex[2])) : qtrk(tpmu_vd[src3+"_qOverP"]->at(src.srcIndex[2]));
	src.q3body  = vtx_charge->at(vtx);
	src.q2body12 = qtrk(trks_qoverp->at(src.trkIndex[0]))+qtrk(trks_qoverp->at(src.trkIndex[1])); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles)
	src.q2body13 = qtrk(trks_qoverp->at(src.trkIndex[0]))+qtrk(trks_qoverp->at(src.trkIndex[2])); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles) 
	src.q2body23 = qtrk(trks_qoverp->at(src.trkIndex[1]))+qtrk(trks_qoverp->at(src.trkIndex[2])); // this is better than the version above since the vertex is built from the id tracks and there's always an id track (also for TrackParticles) 
}