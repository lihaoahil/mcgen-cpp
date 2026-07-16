// Render key mc_analyze histograms to a multi-page PDF for review.
// Usage: root -l -b -q 'validation/plot_hists.C("file.root","out.pdf")'
#include <vector>
#include <string>

void plot_hists(const char* rootfile, const char* pdffile) {
    TFile* f = TFile::Open(rootfile);
    if (!f || f->IsZombie()) {
        printf("cannot open %s\n", rootfile);
        return;
    }
    gStyle->SetOptStat(1110);
    gStyle->SetPalette(kBird);

    // (directory, id, draw option) - one pad each, 4 pads per page.
    struct Item { const char* dir; int id; const char* opt; };
    std::vector<Item> items = {
        {"RAWEVENT", 28, ""},        // beam momentum
        {"RAWEVENT", 8852, ""},      // IM(p2 pbar)
        {"RAWEVENT", 8896, ""},      // -(t - tmin) pbar p
        {"RAWEVENT", 8867, ""},      // pbar cos theta_cm
        {"RAWEVENT", 131, ""},       // p1 lab momentum
        {"RAWEVENT", 132, ""},       // p2 lab momentum
        {"RAWEVENT", 133, ""},       // pbar lab momentum
        {"RAWEVENT", 153, ""},       // CM momentum
        {"RAWEVENT", 8854, "COLZ"},  // Dalitz plot
        {"RAWEVENT", 9870, "COLZ"},  // van Hove
        {"RAWEVENT", 209, "COLZ"},   // proton theta vs p
        {"RAWEVENT", 230, "COLZ"},   // antiproton theta vs p
        {"ACCEPTED", 8852, ""},      // IM(p2 pbar), after acceptance
        {"ACCEPTED", 8896, ""},      // -(t-tmin), after acceptance
        {"ACCEPTED", 8867, ""},      // pbar cos theta_cm, after acceptance
        {"ACCEPTED", 28, ""},        // beam, after acceptance
    };

    TCanvas c("c", "mc_analyze review", 1200, 900);
    c.Divide(2, 2);
    c.Print(Form("%s[", pdffile));
    int pad = 0;
    for (const auto& it : items) {
        TH1* h = nullptr;
        f->GetObject(Form("%s/h%d", it.dir, it.id), h);
        if (!h) {
            printf("missing %s/h%d\n", it.dir, it.id);
            continue;
        }
        c.cd(++pad);
        h->SetTitle(Form("%s: %s", it.dir, h->GetTitle()));
        h->Draw(it.opt);
        if (pad == 4) {
            c.Print(pdffile);
            c.Clear();
            c.Divide(2, 2);
            pad = 0;
        }
    }
    if (pad > 0) c.Print(pdffile);
    c.Print(Form("%s]", pdffile));
    f->Close();
    printf("wrote %s\n", pdffile);
}
