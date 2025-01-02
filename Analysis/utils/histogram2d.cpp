#include <iostream>
#include <map>
#include <utility>
#include <algorithm>
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"

void load_axis_into_vector(const TAxis* axis, std::vector<double>& vector);

class Histogram_2D{
  public:
    Histogram_2D(const char* name, std::vector<double> yaxis, const double xmin, const double xmax);
    Histogram_2D(Histogram_2D& histo) = delete;
    ~Histogram_2D();

    void th2d_add (TH2D*);
    void divide(const Histogram_2D& histo, const std::string& option, const bool);
    void reset();

    void add_x_binning_by_index(const int index, const std::vector<double> xaxis);

    bool can_be_imported(TH2D* histo);
    void print(const char* dir);

    //use this function only for weights, not for counts
    std::shared_ptr<TH2D> get_weights_th2d(const char* name, const char* title);
    std::shared_ptr<TH2D> get_weights_th2d_simpl(const char* name, const char* title);
  private:
    int find_bin_by_value_(const double& y);
    void extrapolate_zero_bins(TH1D* histo);

    std::string name_;

    std::vector<double>                 yaxis_;
    std::vector<std::shared_ptr<TH1D>>  yaxis_content_;
    std::vector<bool> occupancy_;


    double xmin_ = std::numeric_limits<float>::max();
    double xmax_ = std::numeric_limits<float>::lowest();
};

Histogram_2D::~Histogram_2D(){
}

Histogram_2D::Histogram_2D(const char* name, std::vector<double> yaxis, const double xmin, const double xmax){
  yaxis_ = yaxis;
  // std::cout << "start" << std::endl;
  for (std::vector<double>::iterator it = yaxis.begin(); it != std::prev(yaxis.end()); it++){
    // std::cout << *it << std::endl;
    yaxis_content_.push_back(std::make_shared<TH1D>());
    occupancy_.push_back(false);
  }
  // add two more for overflow and underflow bins
  yaxis_content_.push_back(std::make_shared<TH1D>()); occupancy_.push_back(false);
  yaxis_content_.push_back(std::make_shared<TH1D>()); occupancy_.push_back(false);

  xmin_ = xmin;
  xmax_ = xmax;
  name_ = name;
}

void Histogram_2D::add_x_binning_by_index(const int index, const std::vector<double> xaxis){
  std::string name = name_+std::to_string(index);
  if(index > yaxis_content_.size() || index < 0){
    std::runtime_error("Index "+std::to_string(index)+" out of y-axis range");
  }
  if (xaxis.front() != xmax_ || xaxis.back() != xmin_){
    std::runtime_error("Input yaxis min or max values not matching specified min and max values");
  }
  yaxis_content_[index] = std::shared_ptr<TH1D>(new TH1D(name.c_str(), "", xaxis.size()-1, &xaxis[0]));
  occupancy_[index] = true;
}

bool Histogram_2D::can_be_imported(TH2D* histo){
  auto check_axis = [](const std::vector<double>& input_axis, const std::vector<double>& this_axis){
    //if (input_axis.front() != this_axis.front() || input_axis.back() != this_axis.back()) return false;
    bool matching;
    for (auto this_low_edge  : this_axis ){
    for (auto input_low_edge : input_axis){
      matching = (fabs(this_low_edge-input_low_edge) < 2*std::numeric_limits<float>::epsilon());
      if (matching) break;
    } if (!matching) return false;
    }
    return true;
  };

  std::vector<double> input_xaxis;
  std::vector<double> input_yaxis;
  std::vector<double> this_xaxis;

  load_axis_into_vector(histo->GetXaxis(), input_xaxis);
  load_axis_into_vector(histo->GetYaxis(), input_yaxis);

  if (!check_axis(input_yaxis, yaxis_)){
    std::cerr << "Invalid y-axis binning found" << std::endl;
    return false;
  }

  for(int iy = 0; iy < yaxis_.size(); iy++){
    auto this_histo = yaxis_content_[iy];
    load_axis_into_vector(this_histo->GetXaxis(), this_xaxis);
    if(!check_axis(input_xaxis, this_xaxis)){
      std::cerr << "Invalid x axis binning found for y bin n. " << iy << std::endl;
      return false;
    }
    this_xaxis.clear();
  }
  return true;
}

int Histogram_2D::find_bin_by_value_(const double& y){
  for(int iy = 0; iy < yaxis_.size(); iy++){
    if(y < yaxis_[iy] && iy == 0) return iy; // underflow
    if(y >= yaxis_[iy] && iy == (yaxis_.size() - 1)) return iy+1; // overflow
    if(y >= yaxis_[iy] && y < yaxis_[iy+1]) return iy+1;
  }
  throw std::range_error("Value "+std::to_string(y)+" is not in the range of the y axis");
}

void Histogram_2D::th2d_add (TH2D* histo){
  if(!can_be_imported(histo)){
    throw std::invalid_argument("Given TH2D "+std::string(histo->GetName())+" can not be imported");
  }

  auto input_yaxis = histo->GetYaxis();
  auto input_xaxis = histo->GetXaxis();

  for (int iy = 0; iy <= input_yaxis->GetNbins()+1; iy++){
  for (int ix = 0; ix <= input_xaxis->GetNbins()+1; ix++){
    auto bincy = input_yaxis->GetBinCenter(iy);
    auto bincx = input_xaxis->GetBinCenter(ix);
    // std::cout << iy << " " << ix << " " << bincy << " " << bincx << std::endl;
    // if(bincx < xmin_ || bincx >= xmax_ || bincy < yaxis_.front() || bincy >= yaxis_.back()) continue;
    // if(bincy < yaxis_.front() || bincy >= yaxis_.back()) continue;

    auto xhisto = yaxis_content_[find_bin_by_value_(bincy)];
    xhisto->SetBinContent(
      xhisto->FindBin(bincx),
      xhisto->GetBinContent(xhisto->FindBin(bincx)) + histo->GetBinContent(ix, iy));
  }}

  for (auto xhisto : yaxis_content_){
    if (!xhisto->GetSumw2N()){
      xhisto->Sumw2();
    }
  }
}

// void Histogram_2D::extrapolate_zero_bins(TH1D* histo) { // based on density (should be used for the counts)
//     for (int i = 1; i <= histo->GetNbinsX(); ++i) {
//         if (histo->GetBinContent(i) == 0) {
//             double lower_density = 0;
//             double upper_density = 0;
//             double lower_distance = 0;
//             double upper_distance = 0;
//             double lower_error = 0;
//             double upper_error = 0;
//             int count = 0;
//             std::cout << "Extrapolating bin " << i << " " << histo->GetBinContent(i) << std::endl;
//             // Check previous bin
//             if (i > 1) {
//                 double prev_content = histo->GetBinContent(i - 1);
//                 double prev_width = histo->GetBinWidth(i - 1);
//                 if (prev_content != 0) {
//                     lower_density = prev_content / prev_width;
//                     lower_distance = std::abs(histo->GetBinCenter(i) - histo->GetBinCenter(i - 1));
//                     std::cout << "Lower distance: " << lower_distance << std::endl;
//                     std::cout << "Lower density: " << lower_density << std::endl;
//                     lower_error = histo->GetBinError(i - 1) / prev_width;
//                     ++count;
//                 }
//             }

//             // Check next bin
//             if (i < histo->GetNbinsX()) {
//                 double next_content = histo->GetBinContent(i + 1);
//                 double next_width = histo->GetBinWidth(i + 1);
//                 if (next_content != 0) {
//                     upper_density = next_content / next_width;
//                     upper_distance = std::abs(histo->GetBinCenter(i) - histo->GetBinCenter(i + 1));
//                     std::cout << "Upper distance: " << upper_distance << std::endl;
//                     std::cout << "Upper density: " << upper_density << std::endl;
//                     upper_error = histo->GetBinError(i + 1) / next_width;
//                     ++count;
//                 }
//             }

//             if (count > 0) {
//                 double weighted_mean_density = 0;
//                 double weighted_error = 0;
//                 double bin_width = histo->GetBinWidth(i);

//                 if (lower_distance != 0 && upper_distance != 0) {
//                     double total_distance = lower_distance + upper_distance;
//                     weighted_mean_density = (lower_density * upper_distance + upper_density * lower_distance) / total_distance;
//                     weighted_error = std::sqrt(
//                         std::pow((upper_distance / total_distance) * lower_error, 2) +
//                         std::pow((lower_distance / total_distance) * upper_error, 2)
//                     );
//                 } else if (lower_distance != 0) {
//                     weighted_mean_density = lower_density;
//                     weighted_error = lower_error;
//                 } else if (upper_distance != 0) {
//                     weighted_mean_density = upper_density;
//                     weighted_error = upper_error;
//                 }

//                 double real_count = weighted_mean_density * bin_width;
//                 histo->SetBinContent(i, real_count);

//                 // Calculate uncertainty as the propagated error
//                 double uncertainty = weighted_error * bin_width;
//                 histo->SetBinError(i, uncertainty);
//                 std::cout << "Extrapolated bin " << i << " with content " << real_count << " and uncertainty " << uncertainty << std::endl;
//             }
//         }
//     }
// }

void Histogram_2D::extrapolate_zero_bins(TH1D* histo) {
    for (int i = 1; i <= histo->GetNbinsX(); ++i) {
        if (histo->GetBinContent(i) == 0) {
            double lower_value = 0;
            double upper_value = 0;
            double lower_error = 0;
            double upper_error = 0;
            int count = 0;

            // Check previous bin
            if (i > 1) {
                double prev_content = histo->GetBinContent(i - 1);
                if (prev_content != 0) {
                    lower_value = prev_content;
                    lower_error = histo->GetBinError(i - 1);
                    ++count;
                }
            }

            // Check next bin
            if (i < histo->GetNbinsX()) {
                double next_content = histo->GetBinContent(i + 1);
                if (next_content != 0) {
                    upper_value = next_content;
                    upper_error = histo->GetBinError(i + 1);
                    ++count;
                }
            }

            if (count > 0) {
                double mean_value = 0;
                double mean_error = 0;

                if (count == 2) {
                    mean_value = (lower_value + upper_value) / 2.0;
                    mean_error = std::sqrt(std::pow(lower_error, 2) + std::pow(upper_error, 2)) / 2.0;
                } else if (lower_value != 0) {
                    mean_value = lower_value;
                    mean_error = lower_error;
                } else if (upper_value != 0) {
                    mean_value = upper_value;
                    mean_error = upper_error;
                }

                histo->SetBinContent(i, mean_value);
                histo->SetBinError(i, mean_error);
            }
        }
    }
}

void Histogram_2D::divide(const Histogram_2D& histo, const std::string& option, const bool extrapolate){
  auto check_axis = [] (const std::vector<double>& axis1, const std::vector<double>& axis2){
    return axis1 == axis2;
  };

  if (!std::all_of(occupancy_.begin(), occupancy_.end(), [](bool i){return i;})){
    throw std::logic_error("Not all the bins have been initialized");
  }
  if (!check_axis(yaxis_, histo.yaxis_)){
    throw std::logic_error("Invalid x binning detected on denominator for Histogram_2D "+histo.name_);
  }
  
  std::vector<double> thisxaxis;
  std::vector<double> xaxis;
  for (int iy = 0; iy < yaxis_content_.size(); iy++){
    TH1D* thisxhisto = yaxis_content_[iy].get();
    TH1D* xhisto     = histo.yaxis_content_[iy].get();
    load_axis_into_vector(thisxhisto->GetXaxis(), thisxaxis);
    load_axis_into_vector(xhisto->GetXaxis()    , xaxis    );

    if(!check_axis(thisxaxis, xaxis)){
      throw std::logic_error("Invalid x-axis binning found for denominator in y bin n. "+std::to_string(iy)+" for Histogram_2D "+histo.name_);
    }

    if (option == "B"){
      (*thisxhisto).Divide(thisxhisto, xhisto, 1.0, 1.0, "B");
    } else {
      (*thisxhisto).Divide(thisxhisto, xhisto, 1.0, 1.0);
    }

    if (extrapolate && iy != 0 && iy != yaxis_content_.size()-1){
        extrapolate_zero_bins(thisxhisto);
    }
  }
}

std::shared_ptr<TH2D> Histogram_2D::get_weights_th2d(const char* name, const char* title){
  double xwidthmin = std::numeric_limits<float>::max();
  double ywidthmin = std::numeric_limits<float>::max();

  auto get_min_width = [](std::vector<double>& vector, double& min){
    for (int i = 0; i < vector.size()-1; i++) min = std::min(min, vector[i+1] - vector[i]);
  };

  get_min_width(yaxis_, ywidthmin);

  std::vector<double> xaxis;
  for (int iy = 0; iy < yaxis_content_.size(); iy++){
    auto xhisto = yaxis_content_[iy];
    load_axis_into_vector(xhisto->GetXaxis(), xaxis);
    get_min_width(xaxis, xwidthmin);
    xaxis.clear();
  }

  int ny = (yaxis_.back() - yaxis_.front()) / ywidthmin;
  int nx = (xmax_ - xmin_) / xwidthmin;

  auto histo_ = std::make_shared<TH2D>(name, title, nx, xmin_, xmax_, ny, yaxis_.front(), yaxis_.back());

  for(int iy = 1; iy <= ny; iy++){
  for(int ix = 1; ix <= nx; ix++){
    double bincx = histo_->GetXaxis()->GetBinCenter(ix);
    double bincy = histo_->GetYaxis()->GetBinCenter(iy);

    auto xhisto    = yaxis_content_[find_bin_by_value_(bincy)];
    double content = xhisto->GetBinContent(xhisto->FindBin(bincx));
    double error   = xhisto->GetBinError(xhisto->FindBin(bincx));

    histo_->SetBinContent(ix, iy, content);
    histo_->SetBinError  (ix, iy, error  );
  }}

  return histo_;
}

std::shared_ptr<TH2D> Histogram_2D::get_weights_th2d_simpl(const char* name, const char* title){

  std::vector<double> xaxis_;
  auto xhisto = yaxis_content_[0]; // just take first axis and make final histogram in this biining
  load_axis_into_vector(xhisto->GetXaxis(), xaxis_);
  auto histo_ = std::make_shared<TH2D>(name, title, xaxis_.size()-1, &xaxis_[0], yaxis_.size()-1, &yaxis_[0]);
  for(int iy = 0; iy <= histo_->GetNbinsY()+1; iy++){
    for(int ix = 0; ix <= histo_->GetNbinsX()+1; ix++){
      double bincx = histo_->GetXaxis()->GetBinCenter(ix);
      double bincy = histo_->GetYaxis()->GetBinCenter(iy);

      auto xhisto    = yaxis_content_[find_bin_by_value_(bincy)];
      double content = xhisto->GetBinContent(xhisto->FindBin(bincx));
      double error   = xhisto->GetBinError(xhisto->FindBin(bincx));

      histo_->SetBinContent(ix, iy, content);
      histo_->SetBinError  (ix, iy, error  );
  }}
  return histo_;
}

void Histogram_2D::print(const char* dir){
  // TCanvas can;
  for(int iy = 0; iy < yaxis_content_.size(); iy++){
    auto xhisto = yaxis_content_[iy];
    xhisto->Print("all");
    // std::string name = std::string(dir)+"/"+xhisto->GetName()+".pdf";
    // xhisto->Draw("HIST");
    // can.SaveAs(name.c_str(), "pdf");
  }
}

void Histogram_2D::reset(){
  for (int iy = 0; iy < yaxis_content_.size(); iy++){
    auto xhisto = yaxis_content_[iy].get();
    xhisto->Reset();
  }
}

void load_axis_into_vector(const TAxis* axis, std::vector<double>& vector){
  for(int i = 1; i <= axis->GetNbins() + 1; i++) vector.push_back(axis->GetBinLowEdge(i));
}