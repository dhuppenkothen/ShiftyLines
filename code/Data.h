#ifndef _Data_
#define _Data_

#include <vector>

class Data
{
	private:
		std::vector<double> f_left, f_right, y, yerr;
		std::vector<double> f_mid, df;

		const std::vector<double> line_pos;
		const int nlines;

		// Some useful summaries
		double f_min, f_max, f_range;
		void compute_summaries();
	public:
		Data();
		void load(const char* filename);

		// Getters
		const std::vector<double>& get_f_mid() const { return f_mid; }
		const std::vector<double>& get_df() const { return df; }
		const std::vector<double>& get_f_left() const { return f_left; }
		const std::vector<double>& get_f_right() const { return f_right; }
		const std::vector<double>& get_y() const { return y; }
		const std::vector<double>& get_yerr() const { return yerr; }
		const std::vector<double>& get_line_pos() const { return line_pos; }


		double get_f_min() const { return f_min; }
		double get_f_max() const { return f_max; }
		double get_f_range() const { return f_range; }
		int get_nlines() const {return nlines; }
		
	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

#endif

