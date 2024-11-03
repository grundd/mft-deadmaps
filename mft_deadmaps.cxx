#include <ctime>
#include <vector>
#include <map>
#include "CCDB/CcdbApi.h"
//#include "DataFormatsITSMFT/TimeDeadMap.h"

const int n_masked = 28;
int masked_chips[n_masked] = { // as of 2024
  2, 8, 45, 46, 47, 68, 
  102, 103, 104, 161, 162, 163, 
  470, 476, 511, 536, 540, 541, 
  542, 546, 547, 548, 549, 550, 
  551, 579, 580, 581
};
o2::ccdb::CcdbApi api;
long ts_SOR = -1;
long ts_EOR = -1;
long orbit_SOR = -1;
long orbit_EOR = -1;

bool sort_chips (std::tuple<int, float> chip1, std::tuple<int, float> chip2)
{
  if (std::get<1>(chip1) == std::get<1>(chip2)) {
    return std::get<0>(chip1) > std::get<0>(chip2); // sort based on the chip number
  } else {
    return std::get<1>(chip1) > std::get<1>(chip2); // sort based on the "deadness"
  }
}

bool is_masked (int idx_chip)
{
  bool is_masked = false;
  int* masked = std::find(std::begin(masked_chips), std::end(masked_chips), idx_chip);
  if (masked != std::end(masked_chips)) is_masked = true;
  return is_masked;
}

std::string unixts_to_string (long ts_ms, std::string format)
{
  long ts_sec = ts_ms / 1000; // ms to sec
  std::time_t ts_as_timet = ts_sec; // long to time_t
  auto ts_as_tm = std::localtime(&ts_as_timet);
  char buff[80];
  std::strftime(buff, sizeof(buff), format.data(), ts_as_tm);
  std::string date(buff);
  return date;
}

long orbit_to_unixts (long orbit)
{
  if (ts_SOR < 0 || ts_EOR < 0 || orbit_SOR < 0 || orbit_EOR < 0) {
    std::cerr << "Unix timestamps or orbit numbers for SOR+EOR not set correctly!\n";
    return -1;
  }
  // y = a * x + b (y = unixts, x = orbit)
  double a = (double) (ts_EOR - ts_SOR) / (orbit_EOR - orbit_SOR);
  double b = ts_EOR - a * orbit_EOR;
  return (long)(orbit * a + b);
}

template<typename T>
void set_margins (T* c, float t, float r, float b, float l)
{
  c->SetTopMargin(t); 
  c->SetRightMargin(r);
  c->SetBottomMargin(b);
  c->SetLeftMargin(l);
  return;
}

double get_fontsize (int n_chips)
{
  int n1 = 23;
  int n2 = 51;
  double size1 = 0.038;
  double size2 = 0.025;
  double a = (size2 - size1) / (n2 - n1);
  double b = size2 - a * n2;
  return n_chips * a + b;
}

template<typename T>
void draw_axis_graph (TCanvas* c, T* g, std::vector<unsigned long>* orbits)
{
  // to customize:
  int n_labels = 20;

  float y_min = TMath::MinElement(g->GetN(), g->GetY())-1;
  float y_max = TMath::MaxElement(g->GetN(), g->GetY())+1;
  g->GetYaxis()->SetRangeUser(y_min, y_max);
  g->GetYaxis()->SetTickLength(0.01); 
  g->GetYaxis()->SetTitle("#Dead chips");
  // suppress the x-axis labels, ticks and title
  g->GetXaxis()->SetLabelSize(0);
  g->GetXaxis()->SetTitleSize(0);
  g->GetXaxis()->SetTickSize(0);
  double x_min, x_max, y;
  g->GetPoint(0, x_min, y);
  g->GetPoint(g->GetN()-1, x_max, y);
  TAxis *x_ax = g->GetXaxis();
  float edge = (x_max - x_min) * 0.05;
  x_ax->SetLimits(x_min - edge, x_max + edge);

  // x-axis labels and ticks
  c->cd();
  int incr = orbits->size() / (float)(n_labels-1);
  float dy = (y_max - y_min) / (1. - c->GetTopMargin() - c->GetBottomMargin());
  float tick_length = (y_max - y_min) * 0.02;
  for (int i = 0; i < n_labels; i++) 
  {
    double x, y;
    int curr_bin = i * incr;
    g->GetPoint(curr_bin, x, y);
    // tick 
    TLine* l = new TLine(x, y_min, x, y_min + tick_length);
    l->SetLineWidth(2);
    l->Draw();
    // label
    auto orbit = orbits->at(curr_bin);
    std::string label = unixts_to_string(orbit_to_unixts(orbit), "%d/%m, %H:%M");
    float y_pos = y_min - dy * c->GetBottomMargin() / 20;
    TLatex* t = new TLatex(x, y_pos, label.data());
    t->SetTextSize(0.025);
    t->SetTextFont(42);
    t->SetTextAlign(33);
    t->SetTextAngle(45);
    t->Draw();
  }
  return;
}

template<typename T>
void draw_axis_histo (TCanvas* c, T* h, std::vector<unsigned long>* orbits)
{
  // to customize:
  int n_labels = 20;

  h->GetYaxis()->SetTickLength(0.01); 
  h->GetYaxis()->SetLabelSize(get_fontsize(h->GetNbinsY()));
  h->GetYaxis()->SetTitleSize(0);
  // suppress the x-axis labels, ticks and title
  h->GetXaxis()->SetLabelSize(0);
  h->GetXaxis()->SetTitleSize(0);
  h->GetXaxis()->SetTickSize(0);

  // x-axis labels and ticks
  c->cd();
  int incr = orbits->size() / (float)(n_labels-1);
  float y_min = h->GetYaxis()->GetBinLowEdge(1);
  float y_max = h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1);
  float dy = (y_max - y_min) / (1. - c->GetTopMargin() - c->GetBottomMargin());
  float tick_length = (y_max - y_min) * 0.02;
  for (int i = 0; i < n_labels; i++) 
  {
    int curr_bin = i * incr;
    // tick 
    TLine* l = new TLine(curr_bin, y_min, curr_bin, y_min + tick_length);
    l->SetLineWidth(2);
    l->Draw();
    // label
    auto orbit = orbits->at(curr_bin);
    std::string label = unixts_to_string(orbit_to_unixts(orbit), "%d/%m, %H:%M");
    float y_pos = y_min - dy * c->GetBottomMargin() / 20;
    TLatex* t = new TLatex(curr_bin, y_pos, label.data());
    t->SetTextSize(0.025);
    t->SetTextFont(42);
    t->SetTextAlign(33);
    t->SetTextAngle(45);
    t->Draw();
  }
  return;
}

void draw_title (TCanvas* c, std::string title)
{
  c->cd();
  TLatex* t = new TLatex();
  t->SetTextSize(0.032);
  t->SetTextFont(42);
  t->SetTextAlign(22);
  t->DrawLatexNDC(0.5, 0.97, title.data());
  return;
}

void analyze_deadmap (int run, bool verbose = false, bool debug = false)
{
  gStyle->SetOptStat(0);
  int palette[2] = {kWhite, kRed+1};
  gStyle->SetPalette(2, palette);

  auto runInf = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run);
  ts_SOR = runInf.sor;
  ts_EOR = runInf.eor;
  orbit_SOR = runInf.orbitSOR;
  orbit_EOR = runInf.orbitEOR;
  cout << "Aggregated run information:\n"
    << " SOR: " << ts_SOR 
    << " (" << unixts_to_string((long)ts_SOR, "%m/%d/%Y, %H:%M:%S") << "), orbit: " << orbit_SOR << "\n"
    << " EOR: " << ts_EOR 
    << " (" << unixts_to_string((long)ts_EOR, "%m/%d/%Y, %H:%M:%S") << "), orbit: " << orbit_EOR << "\n";

  long duration = (ts_EOR - ts_SOR) / 1e3; // in seconds
  std::cout << " run duration: " << duration / 60 << " min\n";

  std::map<std::string, std::string> metadata;
  metadata["runNumber"] = std::to_string(run);
  auto deadmap = api.retrieveFromTFileAny<o2::itsmft::TimeDeadMap>("MFT/Calib/TimeDeadMap", metadata, ts_SOR);

  std::map<std::string, std::string> headers = api.retrieveHeaders("MFT/Calib/TimeDeadMap", metadata, ts_SOR);
  if (verbose) {
    std::map<std::string, std::string>::iterator it;
    for (it = headers.begin(); it != headers.end(); it++) std::cout << it->first << "\t" << it->second << "\n";
  }

  long val_from = std::stol(headers["Valid-From"]);
  long val_until = std::stol(headers["Valid-Until"]);  

  std::vector<unsigned long> orbits = deadmap->getEvolvingMapKeys();
  unsigned long orbit_first = orbits.front();
  unsigned long orbit_last = orbits.back();
  int n_orbits = orbits.size();

  std::cout << "Deadmap information:\n"
    << " valid from : " << val_from << " (" << unixts_to_string((long)val_from, "%m/%d/%Y, %H:%M:%S") << "\n"
    << " valid until: " << val_until << " (" << unixts_to_string((long)val_until, "%m/%d/%Y, %H:%M:%S") << "\n";
  
  // does the map contain any orbits before SOR and after EOR
  std::cout << " first orbit: " << orbit_first << "\n"
    << "  difference to the orbit at SOR: " << (float)orbit_first-orbit_SOR << " (OK if positive)\n"
    << " last orbit: " << orbit_last << "\n"
    << "  difference to the orbit at EOR: " << (float)orbit_last-orbit_EOR << " (OK if negative)\n"
    << " #points in the deadmap: " << n_orbits << "\n";
  
  int n_chips = 936;
  TGraph* gr_trend_all = new TGraph();
  TGraph* gr_trend_unmasked = new TGraph();
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
    if (orbit > orbit_last) continue; 
    bool consecutive_chips = false;
    uint16_t first_dead = -1;
    int dead_all(0), dead_unmasked(0);

    std::vector<uint16_t> chips = {};
    deadmap->getMapAtOrbit(orbit, chips);
    if (debug) std::cout << "Orbit " << orbit << ", " << chips.size() << " entries:\n";
    for (const auto &chip : chips) 
    {
      if (chip & (uint16_t)(0x8000)) 
      { // bitwise AND: 
        //  only true if "chip" in binary (0b) is 1xxx xxxx xxxx xxxx
        //  the first digit needs to be one, since 0x8000 = 0b1000 0000 0000 0000
        auto first = chip - 0x8000;

        // count the first chip:
        dead_all++;
        if(!is_masked(first)) dead_unmasked++;
        ((TH1C*)arr_chips->At(first))->Fill(i_orb);
        if (debug) std::cout << " " << chip << ":\n" << "  " << first << "\n";

        if (consecutive_chips) std::cerr << "Problem with consecutive chips!\n";
        consecutive_chips = true;
        first_dead = first;
      } else {
        // count all consecutive chips except the last one:
        if (consecutive_chips) { 
          int n_consecutive = chip - first_dead;
          for (int i = 1; i < n_consecutive; i++) {
            int this_chip = first_dead+i;
            dead_all++;
            if(!is_masked(this_chip)) dead_unmasked++;
            ((TH1C*)arr_chips->At(this_chip))->Fill(i_orb);
            if (debug) std::cout << "  " << this_chip << "\n";
          }
        }
        // count the current chip (in the case of consecutive chips, it is the last one):
        dead_all++;
        if(!is_masked(chip)) dead_unmasked++;
        if (debug) std::cout << " " << chip << "\n";
        ((TH1C*)arr_chips->At(chip))->Fill(i_orb);
        consecutive_chips = false;
      }
    }
    if (debug) std::cout << " total dead: " << dead_all << ", dead unmasked: " << dead_unmasked << "\n";
    gr_trend_all->AddPoint(orbit, dead_all);
    gr_trend_unmasked->AddPoint(orbit, dead_unmasked);
    i_orb++;
  }

  // trend of all and unmasked #dead chips
  TCanvas* c0 = new TCanvas("", "", 900, 600);
  set_margins(c0, 0.06, 0.02, 0.13, 0.08);
  gr_trend_all->Draw("AL*");
  draw_axis_graph(c0, gr_trend_all, &orbits);
  draw_title(c0, Form("Run %i: trend of all dead chips", run));
  c0->Print(Form("%i_trend_all.pdf", run));

  TCanvas* c1 = new TCanvas("", "", 900, 600);
  set_margins(c1, 0.06, 0.02, 0.13, 0.08);
  gr_trend_unmasked->Draw("AL*");
  draw_axis_graph(c1, gr_trend_unmasked, &orbits);
  draw_title(c1, Form("Run %i: trend of unmasked dead chips", run));
  c1->Print(Form("%i_trend_unmasked.pdf", run));

  // count chips not present in the deadmap
  int n_ok = 0;
  for (int i_chip = 0; i_chip < n_chips; i_chip++) {
    if (((TH1C*)arr_chips->At(i_chip))->GetEntries() == 0) n_ok++;
    float deadness = ((TH1C*)arr_chips->At(i_chip))->Integral();
    dead_chips.push_back({i_chip, deadness});
  }
  std::cout << "#Chips not present in the dead map: " << n_ok << "\n";

  // sort the chips, the most dead chips first
  std::sort(dead_chips.begin(), dead_chips.end(), sort_chips);
  int n_dead = 0;
  while (std::get<1>(dead_chips[n_dead]) > 0) n_dead++;
  std::cout << "#Chips present in the dead map: " << n_dead << "\n";

  // trend per chip
  int n_unmasked = n_dead - n_masked;
  if (n_unmasked < 0) {
    std:cerr << "Problem with masked chips!\n";
    return;
  }
  TH2C* h_chips_all = new TH2C("", "", n_orbits, 0, n_orbits, n_dead, 0, n_dead);
  TH2C* h_chips_unmasked = new TH2C("", "", n_orbits, 0, n_orbits, n_unmasked, 0, n_unmasked);

  int i_unmasked = 0;
  for (int i_dead = 0; i_dead < n_dead; i_dead++) 
  {
    int idx_chip = std::get<0>(dead_chips[i_dead]);

    // loop over orbits:
    for (int i_orb = 0; i_orb < n_orbits; i_orb++) {
      if (((TH1C*)arr_chips->At(idx_chip))->GetBinContent(i_orb) > 0) {  // the chip was dead in this orbit
        // fill the histogram with all dead chips:
        h_chips_all->Fill(i_orb, n_dead-i_dead-1);
        // fill the histo with unmasked only:
        if (!is_masked(idx_chip)) h_chips_unmasked->Fill(i_orb, n_unmasked-i_unmasked-1);
      } 
    }

    // set custom bin labels
    h_chips_all->GetYaxis()->SetBinLabel(n_dead-i_dead, Form("#%i [%.0f%%]", idx_chip,
      std::get<1>(dead_chips[i_dead]) / n_orbits * 100));
    if (!is_masked(idx_chip)) {
      h_chips_unmasked->GetYaxis()->SetBinLabel(n_unmasked-i_unmasked, Form("#%i [%.0f%%]", idx_chip,
        std::get<1>(dead_chips[i_dead]) / n_orbits * 100));
      i_unmasked++;
    }
  }

  TCanvas* c2 = new TCanvas("", "", 900, 600);
  set_margins(c2, 0.06, 0.02, 0.13, 0.12);
  h_chips_all->Draw("COL");
  draw_axis_histo(c2, h_chips_all, &orbits);
  draw_title(c2, Form("Run %i: all dead chips (total: %i)", run, n_dead));
  c2->Print(Form("%i_dead_all.pdf", run));

  TCanvas* c3 = new TCanvas("", "", 900, 600);
  set_margins(c3, 0.06, 0.02, 0.13, 0.12);
  h_chips_unmasked->Draw("COL");
  draw_axis_histo(c3, h_chips_unmasked, &orbits);
  draw_title(c3, Form("Run %i: unmasked dead chips (total: %i)", run, n_unmasked));
  c3->Print(Form("%i_dead_unmasked.pdf", run));

  return;
}

void mft_deadmaps ()
{
  api.init("http://alice-ccdb.cern.ch");

  analyze_deadmap(559361);

  return;
}