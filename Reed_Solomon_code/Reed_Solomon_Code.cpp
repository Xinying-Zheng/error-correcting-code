#include <bitset>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
using namespace std;

typedef unsigned int uint;

/**
	The main problem is to create a GF(2^6) class contains:
		- add (xor)
		- sub (xor)
		- mul (cycle property)
		- div (similar with mul)

	Then the poly class contains:
		- add (GF.add coeffs of same order terms)
		- sub (GF.sub coeffs of same order terms)
		- mul (pairwise GF.mul coeffs and sum orders)
		- long div (
					go from highest order term of dividend and do GF.div on coeffs
					sub on orders then get p,
					and remainder is dividend - polyMul(p, divisor)
					)

	Relationship of them :
		- poly is a set of GF
		- the ops of poly highly based on GF(coeffs)
*/

class Field {
public:
	/**
		construct 2 looking up tables for mul and div
	*/
	vector<uint> order2coeff;
	vector<uint> coeff2order;
	uint param;

	Field(const uint _param, const uint primitive_poly) : order2coeff(2 * _param), coeff2order(_param), param(_param) {
		/**
			coeff is the 6 bits representation of element of GF(2^6)
			e.g. x^2 = [000100] = 4

			Following procedure will construct the table :
				- order2coeff : 
				given order like x^2 order is 2 and order2coeff[2] = [000100] = 4
				and mul between 2 element may overflow that order at most be 2 * 2^6(param)
				and the ops in GF will mod primitive_poly
			
				- coeff2order :
				this table was created for looking up the order given coeff
				because a element always appears as a 6 bits number
		*/

		uint coeff = 1;
		order2coeff[0] = 1;

		for (uint order = 1; order < 2 * param; ++order) {
			coeff *= 2;
			coeff = coeff >= param ? coeff ^ primitive_poly : coeff;
			order2coeff[order] = coeff;
			
			if (order < param) coeff2order[coeff] = order;
		}
	}

	// add and sub are both XOR
	uint add(const uint l, const uint r) {
		return l ^ r;
	}

	uint sub(const uint l, const uint r) {
		return l ^ r;
	}

	// when mulFromCoeff we given coeffs and a good way to mul
	// is change coeff to order then sum up orders then turn back to coeff
	uint mulFromCoeff(const uint l, const uint r, bool return_coeff=true) {
		// when mul with 0s return 0 since coeff = 0 is invalid
		if (l == 0 || r == 0) return 0;
		uint order = coeff2order[l] + coeff2order[r];
		if (return_coeff) return order2coeff[order];
		return order > 127 ? order - 127 : order;
	}

	// when given orders of 2 terms directly sum up them to check coeff
	uint mulFromOrder(const uint l, const uint r, bool return_coeff = true) {
		uint order = l + r;
		if (return_coeff) return order2coeff[order];
		return order > 63? order - 63 : order;
	}

	// same as mul, change coeff 2 order and sub
	uint divFromCoeff(const uint l, const uint r, bool return_coeff = true) {
		// 0 is invalid order return 0
		if (l == 0 || r == 0) return 0;
		uint order = param - 1 + coeff2order[l] - coeff2order[r];
		if (return_coeff) return order2coeff[order];
		return order > 63 ? order - 63 : order;
	}

	// directly sub orders
	uint divFromOrder(const uint l, const uint r, bool return_coeff = true) {
		uint order = param - 1 + l - r;
		if (return_coeff) return order2coeff[order];
		return order > 63 ? order - 63 : order;
	}

	// this is a easy problem in leetcode 
	// when sum over inputs N times N times XOR is 0 when N is even itself when N is odd
	uint sumOverN(const uint l, const uint N) {
		return N & 1 ? l : 0;
	}

	// info function for debugging
	void info() {
		cout << "order2coeff : " << endl;
		for (auto i : order2coeff	) cout << i <<' ';
		cout << endl;
		cout << "coeff2order : " << endl;
		for (auto i : coeff2order) cout << i <<' ';
		cout << endl;
	}
};

class Poly {
public:
	uint order;
	vector<uint> coeff;
	
	Poly(vector<uint> _coeff) :coeff(_coeff.begin(), _coeff.end()) {
		// check for highest order term
		int highest = 0, len = _coeff.size();
		for (int i = 0; i < len; i++)
			if (_coeff[i] != 0) highest = i;
		order = highest;
	}

	// info function for debugging
	void info() {
		cout << "order is " << order << endl;
		for (uint i = 0; i < coeff.size(); i++) {
			cout << coeff[i] << "x^" << i;
			if (i != coeff.size() - 1) cout << " + ";
		}
		cout << endl;
	}
};

Poly polyAdd(Field& f, const Poly& l, const Poly& r) {
	// merge coeff of same order terms
	uint res_order = l.order > r.order ? l.order : r.order;
	vector<uint> res_coeff(res_order + 1);
	
	for (uint i = 0; i <= res_order; ++i) {
		uint l_coeff = i <= l.order ? l.coeff[i] : 0;
		uint r_coeff = i <= r.order ? r.coeff[i] : 0;
		res_coeff[i] = f.add(l_coeff, r_coeff);
	}

	return Poly(res_coeff);
}

Poly polySub(Field& f, const Poly& l, const Poly& r) {
	// merge coeff of same order terms
	uint res_order = l.order > r.order ? l.order : r.order;
	vector<uint> res_coeff(res_order + 1);

	for (uint i = 0; i <= res_order; ++i) {
		uint l_coeff = i <= l.order ? l.coeff[i] : 0;
		uint r_coeff = i <= r.order ? r.coeff[i] : 0;
		res_coeff[i] = f.sub(l_coeff, r_coeff);
	}

	return Poly(res_coeff);
}

Poly polyMul(Field& f, const Poly& l, const Poly& r) {
	// 2 loops go from low order term to high order term then do f.mul and f.add
	uint res_order = l.order + r.order;
	vector<uint> res_coeff(res_order + 1);

	for (uint i = 0; i <= l.order; ++i) {
		// when mul with 0 return 0
		// if (l.coeff[i] == 0) continue;
		for (uint j = 0; j <= r.order; ++j) {
			// if (r.coeff[j] == 0) continue;
			uint prod_order = i + j;
			uint prod_coeff = f.mulFromCoeff(l.coeff[i], r.coeff[j]);
			res_coeff[prod_order] = f.add(res_coeff[prod_order], prod_coeff);
		}
	}

	return Poly(res_coeff);
}

pair<Poly, Poly> polyMod(Field& f, Poly& dividend, Poly& divisor) {
	// remainder is defined as dividend - polyMul(p, divisor) so initialized same as dividend
	Poly remainder = Poly(dividend.coeff);
	// q at most has the same degree with dividend
	vector<uint> q_coeffs(dividend.order + 1);
	
	// leading order of divisor
	uint divisor_leading = f.coeff2order[divisor.coeff[divisor.order]];
	// long division steps along one order at a time, starting at the highest order
	for (uint i = dividend.order; i > 0; i--) {

		if (i < divisor.order) {
			break;
		}
		if (remainder.coeff[i] == 0) {
			continue;
		}

		// calculate coeff and order of q
		uint q_order = i - divisor.order;
		uint q_coeff = f.divFromOrder(f.coeff2order[remainder.coeff[i]], divisor_leading);
		q_coeffs[q_order] = q_coeff;

		// now that we've chosen q, mul the divisor by q and sub
		for (uint j = 0; j <= divisor.order; j++) {
			if (divisor.coeff[j] == 0) {
				continue;
			}
			// mul and sub
			remainder.coeff[j + q_order] = f.sub(remainder.coeff[j + q_order],
				f.mulFromCoeff(divisor.coeff[j], q_coeff));
		}
	}

	Poly quotient(q_coeffs);

	return { quotient, remainder };
}

Poly polyDer(Field& f, Poly& poly) {
	// if f(x) = a(n)*x^n + ... + a(1)*x + a(0)
	// then f'(x) = n*a(n)*x^(n-1) + ... + 2*a(2)*x + a(1)
	// where n*a(n) = sum(k=1, n, a(n)) e.g. the nth sum of a(n) in GF(2^8)
	// assumes der.order = poly.order - 1
	vector<uint> res_coeff(poly.order);
	for (uint i = 0; i < res_coeff.size(); i++) {
		// we're filling in the ith power of der, so we look ahead one power in poly
		// f(x) = a(i + 1)*x^(i + 1) -> f'(x) = (i + 1)*a(i + 1)*x^i
		// where (i + 1)*a(i + 1) is the sum of a(i + 1) (i + 1) times, not the product
		res_coeff[i] = f.sumOverN(poly.coeff[i + 1], i + 1);
	}

	return Poly(res_coeff);
}

uint polyEval(Field& f, Poly& poly, uint val) {
	// evaluate the polynomial poly at a particular element val
	// a poly looks like a + b*x + c*x^2 ...
	// when evaluate x on a element of field it goes to be sum of mul of elements
	if (val == 0) {
		return poly.coeff[0];
	}

	uint res = 0;

	// we're going to start at 0th order and multiply by val each time
	uint val_exponentiated = f.coeff2order[1];
	uint val_log = f.coeff2order[val];

	for (uint i = 0; i <= poly.order; i++) {
		if (poly.coeff[i] != 0) {
			// multiply-accumulate by the next coeff times the next power of val
			// initially val_exponentiated will be 0 or 63
			// then mulFromOrder take order representation of poly.coeff[i]
			// and 0(or 63) as input then add then return order2coeff[res]
			// since 0 + poly.coeff[i] == 63 + poly.coeff[i](cyclic) so res will
			// only add 0 order term itself so no wrong
			auto temp = f.mulFromOrder(f.coeff2order[poly.coeff[i]], val_exponentiated);
			res = f.add(res, temp);
		}
		// now advance to the next power
		// add ther order together so each time val_exponentiated += val_log
		val_exponentiated = f.mulFromOrder(val_exponentiated, val_log, false);
	}
	return res;
}

Poly encode(Field& f, vector<uint>& info_bits_coeff, vector<uint>& gen_coeff) {
	// step 1 : contruct x_r
	vector<uint> x_r_coeff(22);
	x_r_coeff.back() = 1;
	Poly x_r(x_r_coeff);

	// step 2 : polyMul(f, info_bits, x_r)
	Poly info_bits(info_bits_coeff);
	Poly mul_res = polyMul(f, info_bits, x_r);
	
	// step 3 : construct code
	Poly gen(gen_coeff);
	auto div_res = polyMod(f, mul_res, gen);
	Poly code = polySub(f, mul_res, div_res.second);

	return code;
}

pair<Poly, Poly> polyExgcd(Field& f,Poly& x_r, Poly& S0, uint mu, uint nu) {
	Poly u0({ 1 }), u1({ 0 });
	Poly v0({ 0 }), v1({ 1 });
	Poly r0 = x_r, r1 = S0;

	while (r0.order > nu || v0.order > mu) {
		auto div_res = polyMod(f, r0, r1);
		Poly q = div_res.first;
		// remainder
		Poly temp = div_res.second;
		// change dividend and divisor
		r0 = r1;
		r1 = temp;
		// update u and v
		temp = polySub(f, u0, polyMul(f, q, u1));
		u0 = u1;
		u1 = temp;
		temp = polySub(f, v0, polyMul(f, q, v1));
		v0 = v1;
		v1 = temp;
	}

	return { v0, r0 };
}

string decode(Field& f, vector<string>& code) {
	int code_len = code.size();
	vector<int> stars_indices;
	for (int i = 0; i < code_len; i++) {
		if (code[i] == "*") stars_indices.push_back(i);
	}

	uint e0 = stars_indices.size();
	
	Poly erasure_loc = e0 > 0 ? Poly({ 1, f.order2coeff[stars_indices[0]] }) : Poly({ 1 });
	//cout << "begining erasure_loc info:" << endl;
	//erasure_loc.info();

	for (int i = 1; i < int(e0); i++) {
		Poly sub_term({ 1, f.order2coeff[stars_indices[i]] });
		erasure_loc = polyMul(f, erasure_loc, sub_term);
	}
	//cout << "total erasure_loc info:" << endl;
	//erasure_loc.info();

	vector<uint> R_prim(code.size());
	for (int i = 0; i < code_len; i++) {
		if (code[i] == "*") continue;
		R_prim[i] = stoi(code[i]);
	}

	//for (uint i : R_prim) cout << i << ' ';
	//cout << endl;

	vector<uint> S_coeffs(21);
	for (int j = 1; j < 22; j++) {
		uint sum = 0;
		for (int i = 0; i < 63; i++) {
			int order = (j * i) % 63;
			uint mul_res = f.mulFromCoeff(R_prim[i], f.order2coeff[order]);
			sum = f.add(sum, mul_res);
		}
		S_coeffs[j - 1] = sum;
	}

	Poly S_x(S_coeffs);
	if (S_x.order == 0 && S_x.coeff[0] == 0) {
		cout << "deocode successfully from R_prim" << endl;
		string res;
		for (uint i : R_prim) {
			res += to_string(i) + " ";
		}
		cout << res << endl;
		return res;
	}
	//cout << "S_x info:" << endl;
	//S_x.info();

	vector<uint> x_r_coeff(22);
	x_r_coeff.back() = 1;
	Poly x_r(x_r_coeff);
	//cout << "x_r info : " << endl;
	//x_r.info();

	Poly mul_res = polyMul(f, erasure_loc, S_x);
	auto div_res = polyMod(f, mul_res, x_r);

	Poly S0 = div_res.second;
	//S0.info();

	uint mu = floor((double)(21 - e0) / 2);
	uint nu = ceil((double)(21 + e0) / 2) - 1;
	//cout << "mu is : " << mu << endl;
	//cout << "nu is : " << nu << endl;

	// Poly exgcd
	auto exgcd_res = polyExgcd(f, x_r, S0, mu, nu);
	
	Poly error_loc = exgcd_res.first;
	Poly omega = exgcd_res.second;

	//cout << "error_loc info" << endl;
	//error_loc.info();
	//cout << "omega info" << endl;
	//omega.info();

	// Time-Domain Completion
	vector<uint> E_coeff(63);
	Poly erasure_error_loc = polyMul(f, erasure_loc, error_loc);
	
	//cout << "erasure_error_loc info" << endl;
	//erasure_error_loc.info();

	if (erasure_error_loc.coeff[0] == 0 || omega.order >= e0 + error_loc.order) {
		cout << "can't decode" << endl;
		return "can't decode";
	}
	else {
		int cnt = 0;
		Poly derivative = polyDer(f, erasure_error_loc);
		for (uint i = 0; i < 63; i++){
			uint order = 63 - i;
			uint eval_sigma = polyEval(f, erasure_error_loc, f.order2coeff[order]);
			uint eval_sigma_prim = polyEval(f, derivative, f.order2coeff[order]);
			if (eval_sigma == 0 && eval_sigma_prim != 0) {
				cnt++;
				uint eval_omega = polyEval(f, omega, f.order2coeff[order]);
				uint div_res = f.divFromCoeff(eval_omega, eval_sigma_prim);
				E_coeff[i] = div_res;
			}
			else E_coeff[i] = 0;
		}

		if (cnt == erasure_error_loc.order)
			cout << "decode successfully" << endl;
		else {
			cout << "can't decode" << endl;
			return "can't decode";
		}
	}
	
	Poly E(E_coeff);
	Poly decode_res = polySub(f, R_prim, E);
	string res;
	for (uint i : decode_res.coeff) res += to_string(i) + " ";
	cout << res << endl;
	return res;
}

void split(const string& s, vector<string>& res, const string& delimiter = " ") {
	string::size_type last_pos = s.find_first_not_of(delimiter, 0);
	string::size_type pos = s.find_first_of(delimiter, last_pos);
	while (pos != string::npos || last_pos != string::npos) {
		res.emplace_back(s.substr(last_pos, pos - last_pos));
		last_pos = s.find_first_not_of(delimiter, pos);
		pos = s.find_first_of(delimiter, last_pos);
	}
}

int main() {
	Field f(64, 67);
	
	cout << "Begin decoding ..." << endl;
	string in_file, out_file;
	cout << "Please give input filename : ";
	cin >> in_file;
	cout << "Please give output filename : ";
	cin >> out_file;

	ifstream infs(in_file);
	ofstream outfs(out_file);

	if (!infs.is_open()) {
		cout << "Open file Error !" << endl;
		exit(1);
	}

	string line;

	while (getline(infs, line)){
		cout << line << endl;
		vector<string> codes;
		split(line, codes);

		// decode
		string decode_res = decode(f, codes);
		if (outfs.is_open()) outfs << decode_res + "\n";
		cout << endl;
	}

	infs.close();
	outfs.close();
}