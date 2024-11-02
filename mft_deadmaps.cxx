#include <vector>
#include <map>
#include "CCDB/CcdbApi.h"
//#include "DataFormatsITSMFT/TimeDeadMap.h"

o2::ccdb::CcdbApi api;

bool sort_chips (std::tuple<int, float> chip1, std::tuple<int, float> chip2)
{
  if (std::get<1>(chip1) == std::get<1>(chip2)) {
    return std::get<0>(chip1) > std::get<0>(chip2); // sort based on the chip number
  } else {
    return std::get<1>(chip1) > std::get<1>(chip2); // sort based on the "deadness"
  }
}

long get_timestamp (int run, std::string str) 
{
  long timestamp = -1;
  std::map<std::string, std::string> hdRCT = api.retrieveHeaders("RCT/Info/RunInformation", map<std::string, std::string>(), run);
  const auto startRCT = hdRCT.find(str);
  if (startRCT != hdRCT.end()) {
    timestamp = stol(startRCT->second);
    std::cout << str << " found, timestamp: " << timestamp << "\n";
  } else {
    std::cout << str << " not found in headers!" << "\n";
  }
  return timestamp;
}

bool download_deadmap (int run, bool verbose = true, bool debug = false)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kWaterMelon);

  auto runinf = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run);
  long ts_sor = runinf.sor;
  cout << " SOR: " << ts_sor << "\n";
  //long ts_eor

  /*
  std::cout << "orbitSOR=" << ri.orbitSOR << " orbitEOR=" << ri.orbitEOR << " orbitsPerTF=" << ri.orbitsPerTF  << ri.orbitEOR << " sor=" << ri.sor << " orbitReset=" << ri.orbitReset  << "\n";
  orbitSOR=18001664 orbitEOR=489344256 orbitsPerTF=32489344256 sor=1730309106567 orbitReset=1730307505776379
  ri.grpECS->print();
  */

  long ts_stf = get_timestamp(run, "STF"); // start of the first TF
  long ts_etf = get_timestamp(run, "ETF"); // end of the last TF
  if (verbose) {
    std::cout << "Run " << run << "\n"
      << " STF: " << ts_stf << "\n"
      << " ETF: " << ts_etf << "\n";
  }

  return false;

  std::map<std::string, std::string> metadata;
  metadata["runNumber"] = std::to_string(run);
  auto deadmap = api.retrieveFromTFileAny<o2::itsmft::TimeDeadMap>("MFT/Calib/TimeDeadMap", metadata, ts_stf);

  std::map<std::string, std::string> headers = api.retrieveHeaders("MFT/Calib/TimeDeadMap", metadata, ts_stf);
  if (verbose) {
    std::map<std::string, std::string>::iterator it;
    for (it = headers.begin(); it != headers.end(); it++) std::cout << it->first << "\t" << it->second << "\n";
  }

  long val_from = std::stol(headers["Valid-From"]);
  long val_until = std::stol(headers["Valid-Until"]);  

  std::vector<unsigned long> orbits = deadmap->getEvolvingMapKeys();
  unsigned long first_orbit = orbits.front();
  unsigned long last_orbit = orbits.back();
  std::cout << "First orbit: " << first_orbit << ", last orbit: " << last_orbit << "\n";
  float bin_width = (last_orbit - first_orbit) / 100000;

  long duration = (ts_etf - ts_stf) / 1e3; // duration in secs
  int n_orbits = orbits.size();

  std::cout << "Duration: " << duration / 60 << " mins, #bins: " << n_orbits << "\n"
    << " One bin ~" << (float)duration / n_orbits << " sec\n";
  
  int n_chips = 936;
  TGraph* gr_trend = new TGraph();
  TObjArray* arr_chips = new TObjArray(n_chips);
  arr_chips->SetOwner(true);
  for (int i = 0; i < n_chips; i++)
  {
    TH1C* h_chip = new TH1C(Form("chip_%i", i), "", n_orbits, 0, n_orbits);
    arr_chips->AddAt(h_chip, i);
  }
  std::vector<std::tuple<int, float>> dead_chips;

  int i_orb = 0;
  for (const auto &orbit : orbits) 
  {
    if (orbit > last_orbit) continue; 
    bool consecutive_chips = false;
    uint16_t first_dead = -1;
    int total_dead = 0;

    std::vector<uint16_t> chips = {};
    long lower_orbit = deadmap->getMapAtOrbit(orbit, chips);
    if (debug) std::cout << "Orbit " << orbit << ", " << chips.size() << " entries:\n";
    for (const auto &chip : chips) 
    {
      if (chip & (uint16_t)(0x8000)) 
      { // bitwise and 
        // only true if "chip" in binary (0b) is 1xxx xxxx xxxx xxxx
        // The first digit needs to be one, since 0x8000 = 0b1000 0000 0000 0000
        auto first = chip - 0x8000;
        if (debug) std::cout << " " << chip << ":\n" << "  " << first << "\n";
        ((TH1C*)arr_chips->At(first))->Fill(i_orb);
        if (consecutive_chips) std::cerr << "Problem with consecutive chips!\n";
        consecutive_chips = true;
        first_dead = first;
      } else {
        if (consecutive_chips) {
          int n_consecutive = chip - first_dead + 1;
          total_dead += n_consecutive;
          for (int i = 1; i < n_consecutive-1; i++) 
          {
            if (debug) std::cout << "  " << first_dead+i << "\n";
            ((TH1C*)arr_chips->At(first_dead+i))->Fill(i_orb);
          }
        } else {
          total_dead++;
        }
        if (debug) std::cout << " " << chip << "\n";
        ((TH1C*)arr_chips->At(chip))->Fill(i_orb);
        consecutive_chips = false;
      }
    }
    if (debug) std::cout << " total: " << total_dead << "\n";
    gr_trend->AddPoint(orbit, total_dead);
    i_orb++;
  }

  TCanvas* c1 = new TCanvas("", "", 900, 600);
  gr_trend->Draw("AL");
  int n_ok = 0;
  for (int i_chip = 0; i_chip < n_chips; i_chip++) {
    if (((TH1C*)arr_chips->At(i_chip))->GetEntries() == 0) n_ok++;
    float deadness = ((TH1C*)arr_chips->At(i_chip))->Integral();
    dead_chips.push_back({i_chip, deadness});
  }
  std::cout << "#Chips not present in the dead map: " << n_ok << "\n";

  // sort the chips in decreasing order
  std::sort(dead_chips.begin(), dead_chips.end(), sort_chips);

  int n_dead = 0;
  while (std::get<1>(dead_chips[n_dead]) > 0) n_dead++;
  std::cout << "#Chips present in the dead map: " << n_dead << "\n";

  TH2C* h_trend_chip = new TH2C("", "", n_orbits, 0, n_orbits, n_dead, 0, n_dead);
  for (int i_dead = 0; i_dead < n_dead; i_dead++) {
    int idx_chip = std::get<0>(dead_chips[i_dead]);
    for (int i_orb = 0; i_orb < n_orbits; i_orb++) {
      if (((TH1C*)arr_chips->At(idx_chip))->GetBinContent(i_orb) > 0) h_trend_chip->Fill(i_orb, n_dead-i_dead-1);
    }
    h_trend_chip->GetYaxis()->SetBinLabel(n_dead-i_dead, Form("#%i [%.0f%%]", idx_chip,
      std::get<1>(dead_chips[i_dead]) / n_orbits * 100
    ));
  }

  TCanvas* c2 = new TCanvas("", "", 900, 600);
  h_trend_chip->Draw("COL");
  return false;
}

void mft_deadmaps ()
{
  api.init("http://alice-ccdb.cern.ch");

  download_deadmap(559211);

  return;
}