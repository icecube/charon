// PYTHIA script to generate neutrino production.
// Calculations for interaction/absorption are taken from S.RITZ, D.SECKEL.
// Nuclear Physics B304 (1988) 877-908
// No EW correction added

// Author : Q.R. Liu


#include "Pythia8/Pythia.h"

// Allow assertions
#include <cassert>

using namespace Pythia8;

//==========================================================================

// Pythia generator.
Pythia pythia;
Pythia hadPythia;

// Debug flag.
const bool DEBUG = false;
// Binning definitions.

/* Hist nuAll   ("all flavor spectrum",             nBinIn, 0., xMaxIn); */
Hist nuE;     
Hist nuEBar;  
Hist nuMu;    
Hist nuMuBar;
Hist nuTau;   
Hist nuTauBar;
//==========================================================================

// Pretty print.

template <class T>
void prettyPrint(string name, T value){
  if (!DEBUG) return;
  cout << "||| " << setw(20) << name << " = " << setw(10) << value << endl;
  return;
}

//==========================================================================

// Check if particle is a B Meson.

bool isBMeson(int& id){
  int idAbs = abs(id);
  if (idAbs == 511    || idAbs == 521    || idAbs == 10511  ||
      idAbs == 10521  || idAbs == 513    || idAbs == 523    ||
      idAbs == 10513  || idAbs == 10523  || idAbs == 20513  ||
      idAbs == 20523  || idAbs == 515    || idAbs == 525    ||
      idAbs == 531    || idAbs == 10531  || idAbs == 533    ||
      idAbs == 10533  || idAbs == 20533  || idAbs == 535    ||
      idAbs == 541    || idAbs == 10541  || idAbs == 543    ||
      idAbs == 10543  || idAbs == 20543  || idAbs == 545)
    return true;
  else return false;
}

//==========================================================================

// Check if particle is a B Bbar.

bool isBBBar(int& id){
  int idAbs = abs(id);
  if (idAbs == 551     || idAbs == 10551   || idAbs == 100551 ||
      idAbs == 110551  || idAbs == 200551  || idAbs == 210551 ||
      idAbs == 553     || idAbs == 10553   || idAbs == 20553  ||
      idAbs == 30553   || idAbs == 100553  || idAbs == 110553 ||
      idAbs == 120553  || idAbs == 130553  || idAbs == 200553 ||
      idAbs == 210553  || idAbs == 220553  || idAbs == 300553 ||
      idAbs == 9000533 || idAbs == 9010533 || idAbs == 555    ||
      idAbs == 10555   || idAbs == 20555   || idAbs == 100555 ||
      idAbs == 110555  || idAbs == 120555  || idAbs == 200555 ||
      idAbs == 557     || idAbs == 100557)
    return true;
  else return false;
}

//==========================================================================

// Check if particle is a B Baryon.

bool isBBaryon(int& id){
  int idAbs = abs(id);
  if (idAbs == 5122 || idAbs == 5112 || idAbs == 5212 || idAbs == 5222 ||
      idAbs == 5114 || idAbs == 5214 || idAbs == 5224 || idAbs == 5132 ||
      idAbs == 5232 || idAbs == 5312 || idAbs == 5322 || idAbs == 5314 ||
      idAbs == 5324 || idAbs == 5332 || idAbs == 5334 || idAbs == 5142 ||
      idAbs == 5242 || idAbs == 5412 || idAbs == 5422 || idAbs == 5414 ||
      idAbs == 5424 || idAbs == 5342 || idAbs == 5432 || idAbs == 5434 ||
      idAbs == 5442 || idAbs == 5444 || idAbs == 5512 || idAbs == 5522 ||
      idAbs == 5514 || idAbs == 5524 || idAbs == 5532 || idAbs == 5534 ||
      idAbs == 5542 || idAbs == 5544 || idAbs == 5554)
    return true;
  else return false;
}

//==========================================================================

// Check if particle is a D Meson.

bool isDMeson(int& id){
  int idAbs = abs(id);
  if (idAbs == 411    || idAbs == 421    || idAbs == 431    ||
      idAbs == 10411  || idAbs == 10421  || idAbs == 413    ||
      idAbs == 423    || idAbs == 10413  || idAbs == 10423  ||
      idAbs == 20413  || idAbs == 20423  || idAbs == 415    ||
      idAbs == 425    || idAbs == 10431  || idAbs == 433    ||
      idAbs == 10433  || idAbs == 20433  || idAbs == 435)
    return true;
  else return false;
}

//==========================================================================

// Check if particle is a C Meson.

bool isCBaryon(int& id){
  int idAbs = abs(id);
  if (idAbs == 4123 || idAbs == 4222 || idAbs == 4212 || idAbs == 4112 ||
      idAbs == 4224 || idAbs == 4214 || idAbs == 4114 || idAbs == 4232 ||
      idAbs == 4132 || idAbs == 4322 || idAbs == 4312 || idAbs == 4324 ||
      idAbs == 4314 || idAbs == 4332 || idAbs == 4334 || idAbs == 4412 ||
      idAbs == 4422 || idAbs == 4414 || idAbs == 4424 || idAbs == 4432 ||
      idAbs == 4434 || idAbs == 4444)
    return true;
  else return false;
}

//==========================================================================

// Fill histogram with final state neutrinos.

void fillFSNeutrino(Event& event, double& weight){
  for (int i = 0; i < event.size(); ++i){
    if (event[i].isFinal()) {
      int idFS = event[i].id();
      double eFS = event[i].e();
      if      (idFS ==  14) nuMu.fill(eFS, weight);
      else if (idFS == -14) nuMuBar.fill(eFS, weight);
      else if (idFS ==  12) nuE.fill(eFS, weight);
      else if (idFS == -12) nuEBar.fill(eFS, weight);
      else if (idFS ==  16) nuTau.fill(eFS, weight);
      else if (idFS == -16) nuTauBar.fill(eFS, weight);
    }
  }
  return;
}

//==========================================================================

// Decay a particle.
// Input event record is updated with the decay event record.

bool decay(int& id, Vec4& p4Vec, Event& event){
  if (DEBUG) cout << "||| Entering decay" << endl;

  prettyPrint("id",  id);
  prettyPrint("p4",  p4Vec);

  // Reset event record and set particle to decay.
  pythia.event.reset();
  pythia.particleData.mayDecay(id, true);

  // Decay with Pythia.
  double mass = pythia.particleData.m0(id);
  pythia.event.append(id, 91, 0, 0, p4Vec, mass);
  if (!pythia.moreDecays()) return false;

  // Reset decay of particle.
  pythia.particleData.mayDecay(id, false);

  // Set the event record.
  event = pythia.event;
  if (DEBUG) event.list();

  if (DEBUG) cout << "||| Leaving decay" << endl;
  return true;
}

//==========================================================================

// Hadronize the jet produced from the heavy quark in an interaction.
// Input event record is updated with the interaction event record.

bool interact(int& id, Vec4& p4Vec, Event& event){
  if (DEBUG) cout << "||| Entering interact" << endl;


  // Define interaction kinematics.
  double m  = hadPythia.particleData.m0(id);
  double mc = hadPythia.particleData.m0(4);
  double qE = p4Vec.e();
  if      ( isBMeson  (id) ) qE = qE*0.7;
  else if ( isDMeson  (id) ) qE = qE*0.6*mc/m;
  else if ( isBBaryon (id) ) qE = qE*0.7;
  else if ( isCBaryon (id) ) qE = qE*0.6*mc/m;
  
  if (qE > m){
  double qP = sqrt(qE*qE-m*m);
  prettyPrint("id",  id);
  prettyPrint("p4",  p4Vec);
  prettyPrint("m",   m);
  prettyPrint("qE",  qE);
  prettyPrint("qP",  qP);

  hadPythia.event.reset();
  hadPythia.event.append(id,   91, 0, 0, 0., 0., qP, qE, m);
  hadPythia.next();

  event = hadPythia.event;
  if (DEBUG) event.list();

  if (DEBUG) cout << "||| Leaving interact" << endl;
  return true;
  }
}

//==========================================================================

// Account for interactions as well as decays of a particle.
// Particle must be either a B Meson, D Meson, B Baryon, C Baryon or B Bbar.
// Events are weighted according to the interaction and decay rates.

bool interDecay(int& id, Vec4& p4Vec,string loc="Sun", double inWeight=1.){
  if (DEBUG) cout << "|| Entering interDecay" << endl;
  if (isnan(inWeight)) assert(0);

  // Define decay length, depending on the boost factor.
  double e           = p4Vec.e();
  prettyPrint ("energy", e);
  if (e <= 10.) return true;
  double mass        = pythia.particleData.m0(id);
  double gamma       = e / mass;
  if (gamma<1.) return true;

  double beta        = sqrt(1. - ((1. / gamma) * (1. / gamma)));
  double lifetime    = pythia.particleData.tau0(id);
  double decayLength = lifetime * gamma;

  //Check if particle decays immediately.
  if (decayLength == 0.){
  cout << "Error, particle ID " << id << " should be set to decay" << endl;
   assert (0);
   }

  // Define interaction length, depending on the type of particle.
  // In units of mm/c.
  //double xSec; // σ_26
  //double c           = 2.99792458e11;// mm/s
  double interLength;
  //double xSec; // σ_26
  //if      ( isBMeson  (id) ) interLength = 2.8e-11*c;
  //else if ( isDMeson  (id) ) interLength = 2.8e-11*c;
  //else if ( isBBaryon (id) ) interLength = 1.6e-11*c;
  //else if ( isCBaryon (id) ) interLength = 1.6e-11*c;
  double xSec;	
  if      ( isBMeson  (id) ) xSec = 1.4;
  else if ( isDMeson  (id) ) xSec = 1.4;
  else if ( isBBaryon (id) ) xSec = 2.4;
  else if ( isCBaryon (id) ) xSec = 2.4;
  else {
    cout << "Error, trying to decay " << id << "\n" << endl;
    assert(0);
  }
  double rho ;
  if      (loc == "Sun")    rho = 138.9;     //sun center density 148.9 g/cm3 	
  else if (loc == "Earth")  rho = 13.08849999999999802; //earth center density 13.088 g/cm3
  double N   = 6.0221409e23;
  prettyPrint("rho",  rho);
  prettyPrint("N",  N);
  prettyPrint("xSec",  xSec);
  interLength  = 1./(rho*N*xSec*1e-26*beta)*10.; //mm/c

  double epsilon = 1e-6;
  if (inWeight < epsilon) return true;

  // Weight by the ratio of the interaction and decay rates.
  double decayWeight = inWeight * interLength / (interLength + decayLength);
  double interWeight = inWeight * decayLength / (interLength + decayLength);

  // Do interactions and decays.
  prettyPrint("particle",     id);
  prettyPrint("decayLength",  decayLength);
  prettyPrint("interLength",  interLength);
  prettyPrint("inWeight",     inWeight);
  prettyPrint("decayWeight",  decayWeight);
  prettyPrint("interWeight",  interWeight);
  prettyPrint("energy",       e);
  prettyPrint("gamma",        gamma);
  prettyPrint("beta",         beta);
  prettyPrint("lifetime",     lifetime);
  prettyPrint("decay length", decayLength);

  if (isnan(decayWeight)) assert(0);
  if (isnan(interWeight)) assert(0);

  Event decayEvent, interEvent;

  if ( !decay   (id, p4Vec, decayEvent) ) return false;
  prettyPrint("decayEvent size", decayEvent.size());
  fillFSNeutrino(decayEvent, decayWeight);
  if ( !interact(id, p4Vec, interEvent) ) return false;
  prettyPrint("interEvent size", interEvent.size());
  fillFSNeutrino(interEvent, interWeight);

  // Fill histogram with final state neutrinos.

  // Handle decay products further.
  for (int i = 0; i < decayEvent.size(); ++i){
    if (decayEvent[i].isFinal() && decayWeight!=0.) {
      int idFS = decayEvent[i].id();

      if (isBMeson(idFS) || isBBaryon(idFS) ||
          isDMeson(idFS) || isCBaryon(idFS)){
        Vec4 p4vecFS = decayEvent[i].p();
        if ( !interDecay(idFS, p4vecFS,loc, decayWeight) ) return false;
      }
    }
  }

  // Handle interaction products further.
  for (int i = 0; i < interEvent.size(); ++i){
    if (interEvent[i].isFinal() && interWeight!=0.) {
      int idFS = interEvent[i].id();

      if (isBMeson(idFS) || isBBaryon(idFS) ||
          isDMeson(idFS) || isCBaryon(idFS)){
        Vec4 p4vecFS = interEvent[i].p();
        if ( !interDecay(idFS, p4vecFS, loc,interWeight) ) return false;
      }
    }
  }

  if (DEBUG) cout << "|| Leaving interDecay" << endl;
  return true;
}

//==========================================================================

// A derived class for (e+ e- ->) GenericResonance -> various final states.
class Sigma1GenRes : public Sigma1Process {

  public:

    // Constructor.
    Sigma1GenRes() {}

    // Evaluate sigmaHat(sHat): dummy unit cross section.
    virtual double sigmaHat() {return 1.;}

    // Select flavour. No colour or anticolour.
    virtual void setIdColAcol() {setId( -11, 11, 999999);
      setColAcol( 0, 0, 0, 0, 0, 0);}

    // Info on the subprocess.
    virtual string name()    const {return "GenericResonance";}
    virtual int    code()    const {return 9001;}
    virtual string inFlux()  const {return "ffbarSame";}

};

//==========================================================================

// Main program function.
  
int main(int argc, char** argv) {

  // A class to generate the fictitious resonance initial state.
  SigmaProcess* sigma1GenRes = new Sigma1GenRes();

  // Hand pointer to Pythia.
  pythia.setSigmaPtr(sigma1GenRes);

  string channel      = string(argv[1]);
  double xMaxIn       = atof(argv[2]); 
  int    nBinIn       = atoi(argv[3]);
  int    seed         = atoi(argv[4]);
  double bin          = double(nBinIn); 
  string location     = string(argv[5]);

  nuE     = Hist("electron neutrino spectrum",      nBinIn, 0., xMaxIn);
  nuEBar  = Hist("anti electron neutrino spectrum", nBinIn, 0., xMaxIn);
  nuMu    = Hist("muon neutrino spectrum",          nBinIn, 0., xMaxIn);
  nuMuBar = Hist("anti muon neutrino spectrum",     nBinIn, 0., xMaxIn);
  nuTau   = Hist("tau neutrino spectrum",           nBinIn, 0., xMaxIn);
  nuTauBar= Hist("anti tau neutrino spectrum",      nBinIn, 0., xMaxIn);
  
  // Read in the rest of the settings and data from a separate file.
  
  //pythia.readFile("channel.cmnd");
  pythia.readFile(channel+"_"+std::to_string(int(xMaxIn))+".cmnd");
  pythia.readFile("decays.cmnd");
  hadPythia.readFile("decays.cmnd");
  // Intialize hadronization Pythia object.
  hadPythia.readString("ProcessLevel:all = off");
  hadPythia.readString("PDF:lepton = off");
  //hadPythia.readString("SoftQCD:all = on");
  hadPythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  hadPythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = "+std::to_string(seed));
  hadPythia.readString("Random:seed = "+std::to_string(seed));
  //hadPythia.readString("WeakBosonExchange:all = on");
  //hadPythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  //hadPythia.readString("WeakSingleBoson:all = on");
  //hadPythia.readString("WeakDoubleBoson:all = on");
  //hadPythia.readString("WeakBosonAndParton:all = on");

  // Initialization.
  if (!DEBUG) pythia.readString("Print:quiet = on");
  if (!DEBUG) hadPythia.readString("Print:quiet = on");
  pythia.init();
  hadPythia.init();

  // Extract settings to be used in the main program.
  int nEvent  = pythia.mode("Main:numberOfEvents");
  int nAbort  = pythia.mode("Main:timesAllowErrors");

  // Progress bar settings.
  int barWidth = 70;

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Progress bar
    cout << "[";
    double progress = (iEvent / (double)nEvent);
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();

    // Generate events. Quit if many failures.
    if (DEBUG) cout << "======================================" << endl;
    if (DEBUG) cout << "GENERATING NEW EVENT #" << iEvent+1 << endl;
    if (!pythia.next()) {
    if (++iAbort < nAbort){ 
	    --iEvent;
	    continue;
    }
      //cout << " Event generation aborted prematurely, owing to error!\n";
      //assert(0);
    }

    // Get the Event Record.
    Event event = pythia.event;
    if (DEBUG) cout << "Event size = " << event.size() << endl;
    if (DEBUG) cout << "======================================" << endl;

    // Loop over all particles and analyze the final-state ones.
    for (int i = 0; i < event.size(); ++i){
      if (event[i].isFinal()) {
        int idFS   = event[i].id();
        double eFS = event[i].e();
        // Fill final state neutrinos into histogram.
        if      (idFS ==  14) nuMu.fill(eFS);
        else if (idFS == -14) nuMuBar.fill(eFS);
        else if (idFS ==  12) nuE.fill(eFS);
        else if (idFS == -12) nuEBar.fill(eFS);
        else if (idFS ==  16) nuTau.fill(eFS);
        else if (idFS == -16) nuTauBar.fill(eFS);

        // Undergo further interactions for particles which interact.
        if (isBMeson(idFS) || isBBaryon(idFS) || isDMeson(idFS) ||
            isCBaryon(idFS)){
          Vec4 p4Vec = event[i].p();
          if ( !interDecay(idFS, p4Vec,location) ){
            if (++iAbort < nAbort) continue;
            cout << " Event generation aborted prematurely, owing to error!\n";
            assert(0);
          };
        }
      }
    }
    // End of event loop.
  }

  // Final statistics and histograms.
  // Rescale according to number of events and binning.
  pythia.stat();
  
  nuE      *= 1. / (nEvent * xMaxIn / bin);
  nuEBar   *= 1. / (nEvent * xMaxIn / bin);
  nuMu     *= 1. / (nEvent * xMaxIn / bin);
  nuMuBar  *= 1. / (nEvent * xMaxIn / bin);
  nuTau    *= 1. / (nEvent * xMaxIn / bin);
  nuTauBar *= 1. / (nEvent * xMaxIn / bin);
  cout << nuMu << nuMuBar << nuE << nuEBar << nuTau << nuTauBar << endl;

  HistPlot hpl("plot_"+channel+"_"+std::to_string(int(xMaxIn))+"_"+location);
  hpl.frame( channel+"_"+std::to_string(int(xMaxIn))+"_"+std::to_string(seed)+"_"+location, "Particle energy spectra", "$E$ (GeV)",
      "$\\mathrm{d}N / \\mathrm{d}E$ (GeV$^{-1}$)");
  hpl.add(nuE,      "-", "$\\nu_e$");
  hpl.add(nuEBar,   "-", "$\\bar{\\nu}_e$");
  hpl.add(nuMu,     "-", "$\\nu_\\mu$");
  hpl.add(nuMuBar,  "-", "$\\bar{\\nu}_\\mu$");
  hpl.add(nuTau,    "-", "$\\nu_\\tau$");
  hpl.add(nuTauBar, "-", "$\\bar{\\nu}_\\tau}$");
  hpl.plot(false); 

  delete sigma1GenRes;
  return 0;
}

//==========================================================================
