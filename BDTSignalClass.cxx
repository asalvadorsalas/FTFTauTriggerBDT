#define BDTSignalClass_cxx
#include "BDTSignalClass.h"
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"

BDTSignalClass::BDTSignalClass(const char* inputfile, const char* outputfile) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TFile *f = new TFile(inputfile,"UPDATE");
  
  f->GetObject("tree_vec",tree);
  
  Init(tree);
}

BDTSignalClass::~BDTSignalClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}
void BDTSignalClass::DeclareHistos(){};
void BDTSignalClass::Loop()
{
if (fChain == 0) return;
//counters 
Int_t Tau12IMcounter =0;
Int_t taucounter =0;
Int_t tautotalcounter = 0;
Int_t missmatchcounter = 0;
Int_t offltrack_total = 0;
Int_t offltrack_lost = 0;
Int_t taucounter_1p = 0;
Int_t taucounter_3p = 0;
Int_t taucounter_1p0n = 0;
Int_t taucounter_1p1n = 0;
Int_t taucounter_1prn = 0;
Int_t offltrack_lost_1p = 0;
Int_t offltrack_lost_3p = 0;
Int_t offltrack_lost_1p0n = 0;
Int_t offltrack_lost_1p1n = 0;
Int_t offltrack_lost_1prn = 0;

//signal track definitions
std::vector<int>*   vCoretrack_signalv0 = 0;
std::vector<int>*   vCoretrack_signalv1 = 0;
TBranch* newBranch0 = tree->Branch("Coretrack_signalv0",&vCoretrack_signalv0);
TBranch* newBranch1 = tree->Branch("Coretrack_signalv1",&vCoretrack_signalv1);

// Additional variables, global
std::vector<float>*   vCoretrack_CaloHadpt = 0;
std::vector<float>*   vCoretrack_CaloEMpt = 0;
std::vector<float>*   vCoretrack_Nprong = 0;
std::vector<float>*   vCoretrack_Nneutral = 0;
std::vector<float>*   vCoretrack_Noffltrk = 0;
std::vector<float>*   vCoretrack_taucharge = 0;
std::vector<float>*   vCoretrack_CaloptCalib = 0;
TBranch* b_tau_ptHad = tree->Branch("Coretrack_CaloHadpt",&vCoretrack_CaloHadpt);
TBranch* b_tau_ptEM = tree->Branch("Coretrack_CaloEMpt",&vCoretrack_CaloEMpt);
TBranch* b_tau_Nprong = tree->Branch("Coretrack_Nprong",&vCoretrack_Nprong);
TBranch* b_tau_Nneutral = tree->Branch("Coretrack_Nneutral",&vCoretrack_Nneutral);
TBranch* b_tau_Noffltrk = tree->Branch("Coretrack_Noffltrk",&vCoretrack_Noffltrk);
TBranch* b_tau_taucharge = tree->Branch("Coretrack_taucharge",&vCoretrack_taucharge);
TBranch* b_tau_CaloptCalib = tree->Branch("Coretrack_CaloptCalib",&vCoretrack_CaloptCalib);

Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes = 0, nb = 0;
int barWidth = 20;
std::cout << "Adding branches " << nentries <<std::endl;
for (Long64_t jentry=0; jentry<nentries;jentry++) { //progress bar
	float progress = float(jentry)/float(nentries);
  if (int(progress*100.0) % 1 ==0){
  	int pos = int(barWidth* progress); std::cout<<"[";
    for (int i=0; i < barWidth; i++){
	 		if(i<pos) std::cout <<"=";
	 		else if (i==pos) std::cout<<">";
	 		else std::cout<<" ";
  	}
      std::cout << "] " << int(progress * 100.0) << "%\r";
      std::cout.flush();
  }
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;

  // Clean variables at the start of the loop
  vCoretrack_signalv0->clear();
  vCoretrack_signalv1->clear();
  vCoretrack_CaloHadpt->clear();
  vCoretrack_CaloEMpt->clear();
  vCoretrack_Nprong->clear();
  vCoretrack_Nneutral->clear();
  vCoretrack_Noffltrk->clear();
  vCoretrack_taucharge->clear();
  vCoretrack_CaloptCalib->clear();
  OfflineTaus.clear();

  for(auto index: *Coretrack_pt ){
    vCoretrack_signalv0->push_back(-1); //init signal flag to -1 to not be used in training
    vCoretrack_signalv1->push_back(-1);
    vCoretrack_CaloHadpt->push_back(-1.);
    vCoretrack_CaloEMpt->push_back(-1.);
    vCoretrack_Nprong->push_back(-1);
    vCoretrack_Nneutral->push_back(-1);
    vCoretrack_Noffltrk->push_back(-1);
    vCoretrack_taucharge->push_back(-99);
    vCoretrack_CaloptCalib->push_back(-1);
  }
  //std::cout << EventNumber << ": " << Coretrack_pt->size() << "/" << vCoretrack_signalv1->size() << std::endl;
  if(passL1_Tau12IM>0){
    Tau12IMcounter++;
    for( int iTau = 0; iTau < tau_offl_pt->size(); ++iTau ) {
      if (tau_offl_isMediumRNN->at(iTau)){ //we only want ismediumRNN

        TLorentzVector tmpTau;
        tmpTau.SetPtEtaPhiM(tau_offl_pt->at(iTau)*0.001,
          tau_offl_eta->at(iTau),
          tau_offl_phi->at(iTau),
          0);  

        ClassTau newTau( tmpTau );
        newTau.Nprong=tau_truth_Nprong->at(iTau);
        newTau.Nneutral=tau_truth_Nneutral->at(iTau);
        newTau.Ntracks=tau_nTracks->at(iTau);  //from preselected HLT
        newTau.HLT_Pt = tau_pt->at(iTau);
        newTau.Calo_Pt = tau_ptTrigCalo->at(iTau)*0.001;
        newTau.Calo_Eta = tau_etaTrigCalo->at(iTau);
        newTau.Calo_Phi = tau_phiTrigCalo->at(iTau);
        newTau.Calo_Had = tau_ptHad->at(iTau);
        newTau.Calo_EM = tau_ptEM->at(iTau);
        newTau.charge=tau_offl_charge->at(iTau);
        newTau.PtCalib=tau_offl_ptFinalCalib->at(iTau);

        std::vector<ClassTrack> tmptracks; //for each tau, save the ass track
        for( int iTrk=0; iTrk < offltrack_pt->at(iTau).size(); ++iTrk){

          TLorentzVector tmpTrk;
          tmpTrk.SetPtEtaPhiM(offltrack_pt->at(iTau).at(iTrk)*0.001,
                              offltrack_eta->at(iTau).at(iTrk),
                              offltrack_phi->at(iTau).at(iTrk),
                              0);  
          ClassTrack newTrk( tmpTrk );
          newTrk.z0=offltrack_z0->at(iTau).at(iTrk);
          newTrk.d0=offltrack_d0->at(iTau).at(iTrk);

          newTrk.nPiHits=offltrack_nPiHits->at(iTau).at(iTrk);
          newTrk.nSiHoles=offltrack_nSiHoles->at(iTau).at(iTrk);  
          newTrk.charge=offltrack_charge->at(iTau).at(iTrk);
          if(offltrack_truthtype->at(iTau).at(iTrk)==1){
            //std::cout << "Tau, trk: " << iTau << " " << iTrk << " "<< newTrk.Pt << std::endl;
            newTrk.istruth=true;
            newTau.hastruthtrk=true;
          }
          tmptracks.push_back(newTrk);
        }

        sort ( tmptracks.begin(), tmptracks.end(), ComparePt() );
        newTau.offlTracks=tmptracks;
        OfflineTaus.push_back(newTau);
      }
    }
    sort( OfflineTaus.begin(), OfflineTaus.end(), ComparePt() );

    for( auto iTau: OfflineTaus ) {  //loop per tau
      tautotalcounter++;
      if( Coretrack_pt->size() > 0 && OfflineTaus.size() > 0 && iTau.Calo_Pt >=0 && iTau.Nprong !=-1 && iTau.offlTracks.size()> 0 && iTau.hastruthtrk ){
        taucounter++;
        if (iTau.Nprong==1){
          taucounter_1p++;
          if(iTau.Nneutral==0) taucounter_1p0n++;
            else if (iTau.Nneutral==1) taucounter_1p1n++;
            else if (iTau.Nneutral>1) taucounter_1prn++;
        }
        else if(iTau.Nprong==3) taucounter_3p++;
      
        //find offl track to match the signal track
        TLorentzVector TCalo; //set HLTCalo
        TCalo.SetPtEtaPhiM(iTau.Calo_Pt, iTau.Calo_Eta, iTau.Calo_Phi, 0);
        ClassTrack targetOfflTrk;

        Bool_t OfflTrkfound = false;
        for ( auto itrk: iTau.offlTracks ) {
          if ( itrk.istruth && fabs(TCalo.Eta()-itrk.Eta) < 0.1 && fabs(TCalo.DeltaPhi(itrk.FourVector)) < 0.1 && itrk.z0 < 225.){
            offltrack_lost++;
            if (iTau.Nprong==1){
              offltrack_lost_1p++;
              if(iTau.Nneutral==0) offltrack_lost_1p0n++;
              else if (iTau.Nneutral==1) offltrack_lost_1p1n++;
              else if (iTau.Nneutral>1) offltrack_lost_1prn++;
            }
            else if(iTau.Nprong==3) offltrack_lost_3p++;
          }
          if( itrk.istruth ){
            targetOfflTrk = itrk;
            OfflTrkfound = true;
            break;
          }
        }
        if(!OfflTrkfound){
          std::cout << "tau w/o truth track, should not happen" << std::endl;
          break;
        } 

        Int_t ntruthtrks = 0;
        for ( auto itrk: iTau.offlTracks ) {
          if (itrk.istruth) ntruthtrks++;
        }

        float tminDR = 9999; //v1,2,
        float tminDz = 9999; //v7
        Int_t stg1index = -1; //v1,v2
        Int_t stg1index001 = -1; // v
        Int_t stg1indexdzdr = -1;
        Int_t stg1indexDz = -1;
        float tminDzdr = 9999; //v6
        ClassTrack target1stg_minDR;
        for( unsigned int iTrk = 0; iTrk < Coretrack_pt->size(); ++iTrk ){

          TLorentzVector tmp1stgtrk;
          tmp1stgtrk.SetPtEtaPhiM(Coretrack_pt->at(iTrk)*0.001,
          Coretrack_eta->at(iTrk),
          Coretrack_phi->at(iTrk),
          0);  
          ClassTrack i1stgtrk( tmp1stgtrk);
          i1stgtrk.z0=Coretrack_z0->at(iTrk);

          if ( fabs(TCalo.Eta() - i1stgtrk.Eta) > 0.1 || fabs(TCalo.DeltaPhi(i1stgtrk.FourVector)) > 0.1 || fabs(i1stgtrk.z0) > 225. ) continue;
          
          vCoretrack_signalv0->at(iTrk)=0;
          vCoretrack_signalv1->at(iTrk)=0;
  
          vCoretrack_CaloHadpt->at(iTrk)=iTau.Calo_Had;
          vCoretrack_CaloEMpt->at(iTrk)=iTau.Calo_EM;          
          vCoretrack_Nprong->at(iTrk)=iTau.Nprong;
          vCoretrack_Nneutral->at(iTrk)=iTau.Nneutral;
          vCoretrack_Noffltrk->at(iTrk)=iTau.offlTracks.size();
          vCoretrack_taucharge->at(iTrk)=iTau.charge;
          vCoretrack_CaloptCalib->at(iTrk)=iTau.PtCalib;

          float DR = targetOfflTrk.FourVector.DeltaR(i1stgtrk.FourVector);
          float Dz = fabs(targetOfflTrk.z0-i1stgtrk.z0);
          if ( DR < 0.05 ){
            vCoretrack_signalv0->at(iTrk)=-1; //veto 0.05
            if( DR < tminDR ){ // minDR
              tminDR = DR;
              target1stg_minDR = i1stgtrk;
              stg1index=iTrk;
            }
            if( Dz < tminDzdr ){ //minDZ within DR 0.05 cone
              tminDzdr = Dz;
              stg1indexdzdr=iTrk;
            }
          }
          if( Dz < 7 && Dz < tminDz){ // minDZ within bound
            tminDz=Dz;
            stg1indexDz=iTrk;
          }
        }//1stg loop

        if (stg1index!=-1){
          vCoretrack_signalv0->at(stg1index)=1; //minDR 0.05 (veto)
          vCoretrack_signalv1->at(stg1index)=1; //minDR 0.05 (only trks inside Calo window) + later 3p exception
        }
      
        //3prong exception for v1
        if (iTau.Nprong==3 && ntruthtrks>1){ 
            std::vector<int> stg1trk_3pindexes;
            for (auto itrk: iTau.offlTracks ){
            if (itrk.istruth){
                stg1trk_3pindexes.push_back(-1);
                //inints
                float minDR = 9999;
                for(unsigned int iTrk =0; iTrk < Coretrack_pt->size(); ++iTrk){
                TLorentzVector tmp1stgtrk;
                tmp1stgtrk.SetPtEtaPhiM(Coretrack_pt->at(iTrk)*0.001,
                Coretrack_eta->at(iTrk),
                Coretrack_phi->at(iTrk),
                0);  
                ClassTrack i1stgtrk( tmp1stgtrk);
                i1stgtrk.z0=Coretrack_z0->at(iTrk);

                if ( fabs(TCalo.Eta() - i1stgtrk.Eta) > 0.1 || fabs(TCalo.DeltaPhi(i1stgtrk.FourVector)) > 0.1 || fabs(i1stgtrk.z0) > 225. ) continue;
            
                float DR = itrk.FourVector.DeltaR(i1stgtrk.FourVector);
                if ( DR < 0.05 ){
                    if ( DR < minDR){
                    minDR = DR;
                    stg1trk_3pindexes.back() = iTrk;
                    }
                }
                }//coretrack loop
            }//istruth requirement
            }//offl trk loop
            if (stg1trk_3pindexes.size()>0 && stg1trk_3pindexes.at(0)!=-1 && vCoretrack_signalv1->at(stg1trk_3pindexes.at(0))!=1)
                std::cout << "ADRI: Dont understand this" << std::endl;
            for(auto indexset: stg1trk_3pindexes)
            if (indexset != -1) vCoretrack_signalv1->at(indexset)=1;

            /*Tried to resolve the overlapping 
            for( unsigned int i=0; i < stg1trk_3pindexes.size(); i++){
                for (unsigned int j=0; j < stg1trk_3pindexes.size(); j++){
                    if (j<=i || stg1trk_3pindexes.at(i)==-1 || stg1trk_3pindexes.at(j)==-1) {
                    continue;
                    }
                    
                    if(stg1trk_3pindexes.at(i)==stg1trk_3pindexes.at(j)){
                    std::cout << "Something wrong: offltrk "<< i << " index " << stg1trk_3pindexes.at(i) << std::endl;
                    std::cout << stg1trk_mindR.at(i) << " vs " << stg1trk_mindR.at(j) << std::endl;
                    int itocheck = i; 
                    if (stg1trk_mindR.at(i) < stg1trk_mindR.at(j)){
                        itocheck = j;
                    }
                    std::vector<float> v = offltrk_mindRs.at(itocheck);
                    sort(v.begin(), v.end());
                    for(unsigned int k=0; k < v.size(); k++)
                        if (v.at(k)!=-1) std::cout << "Init " << k << " " << v.at(k) << std::endl;
                    for(unsigned int k=1; k < v.size() && v.at(k)!=-1; k++) 
                        std::cout  << " sub:" << v.at(k) << std::endl;
                        
                    }
                }
            } */
        } //3prong exception 
      }//tau selection
    }//tau loop
  }
  newBranch0->Fill();
  newBranch1->Fill();

  b_tau_ptHad->Fill();
  b_tau_ptEM->Fill();
  b_tau_Noffltrk->Fill();
  b_tau_taucharge->Fill();
  b_tau_CaloptCalib->Fill();
  b_tau_Nprong->Fill();
  b_tau_Nneutral->Fill();

}//end entryloop  
  newBranch0->Print();
  newBranch1->Print();
 
  b_tau_ptHad->Print();
  b_tau_ptEM->Print();
  b_tau_Noffltrk->Print();
  b_tau_taucharge->Print();
  b_tau_CaloptCalib->Print();
  b_tau_Nprong->Print();
  b_tau_Nneutral->Print();

  tree->Write();

  std::cout << "Total of " << Tau12IMcounter << "events passing L1_Tau12IM" << std::endl;
  std::cout << "Total of " << taucounter << "/"<< tautotalcounter << " tau offl taus considered for training" << std::endl;
  std::cout << "--> split in: 1p,1p0n,1p1n,1pXn,3prong " << taucounter_1p << " " << taucounter_1p0n << " " << taucounter_1p1n << " " << taucounter_1prn << " " << taucounter_3p  << std::endl;
  std::cout << "Total of " << offltrack_lost << " offltrks inside HLTCalo used" << std::endl;
  std::cout << "--> split in: 1p,1p0n,1p1n,1pXn,3prong " << offltrack_lost_1p << " ("<< 100.*(float(taucounter_1p)-float(offltrack_lost_1p))/float(taucounter_1p)  << "%) " << offltrack_lost_1p0n << " ("<< 100.*(float(taucounter_1p0n)-float(offltrack_lost_1p0n))/float(taucounter_1p0n)  << "%) " << offltrack_lost_1p1n << " ("<< 100.*(float(taucounter_1p1n)-float(offltrack_lost_1p1n))/float(taucounter_1p1n)  << "%) " << offltrack_lost_1prn << " ("<< 100.*(float(taucounter_1prn)-float(offltrack_lost_1prn))/float(taucounter_1prn)  << "%) " << offltrack_lost_3p << " ("<< 100.*(float(taucounter_3p)-float(offltrack_lost_3p))/float(taucounter_3p)  << "%) " << std::endl;
}

