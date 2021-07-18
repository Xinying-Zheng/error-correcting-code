#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <list>
#include <bitset>

using namespace std;
// for concise notation
typedef unsigned long long ULL;

class AWGN {
	/**
		An AWGN variable generator follow the instruction
		param:
			SEED : random seed
			RANV : random value
			RANI : random value
	*/
private:
	ULL SEED;
	ULL RANV;
	int RANI;

public:
	AWGN(ULL _SEED) : SEED(_SEED), RANV(_SEED), RANI(0) {}

	double Ranq1() {
		if (RANI == 0) {
			RANV = SEED ^ 4101842887655102017LL;
			RANV ^= RANV >> 21;
			RANV ^= RANV << 35;
			RANV ^= RANV >> 4;
			RANV = RANV * 2685821657736338717LL;
			RANI++;
		}
		RANV ^= RANV >> 21;
		RANV ^= RANV << 35;
		RANV ^= RANV >> 4;
		return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
	}

	pair<double, double> normal(double sigma) {
		double x1, x2, s;

		do {
			x1 = Ranq1();
			x2 = Ranq1();
			//cout << x1 << ' ' << x2 << endl;
			x1 = 2 * x1 - 1;
			x2 = 2 * x2 - 1;
			s = x1 * x1 + x2 * x2;
		} while (s >= 1.0);

		double factor = sqrt(-2 * log(s) / s);

		// lambda function : sigma * x * factor
		auto getRes = [&](double in) {
			return sigma * in * factor;
		};

		return { getRes(x1), getRes(x2) };
	}

	void printInfo() {
		cout << SEED << ' ' << RANV << endl;
	}
};


class Encoder {
	/**
		An encoder generates N + 31 2bits output which comes from
		a generator matrix and AWGN channel

		param:
			info_bits : 63 periodic bits array
			states : 6 bits array
			noise : random variable generator
			N : info_bits length
			ptr : for recycle the info_bits
	*/
private:
	vector<int> info_bits;
	vector<int> states;
	AWGN noise;
	int N, ptr;

public:
	Encoder(int _N, ULL _SEED) : noise(_SEED), N(_N), ptr(0) {
		// info_bits initially set to be 100000
		info_bits = vector<int>(6, 0);
		info_bits[0] = 1;
		// states initially set to be 000000
		states = vector<int>(6, 0);

		// generate info_bits in one period
		while (info_bits.size() < 63) genOneStep();
	}

	// generate the output of generator matrix
	pair<int, int> generate(int bit) {
		// generator matrix operations
		int x1 = bit ^ states[1] ^ states[2] ^ states[4] ^ states[5];
		int x2 = bit ^ states[0] ^ states[1] ^ states[2] ^ states[5];

		// shift left and add a new head
		for (int i = states.size() - 1; i > 0; i--) states[i] = states[i - 1];
		states[0] = bit;

		// map 0 -> 1 and 1 -> -1
		auto mapping = [](int x) { return x == 0 ? 1 : -1; };
		auto res = make_pair(mapping(x1), mapping(x2));

		return res;
	}

	pair<double, double> passAWGN(pair<int, int>& in, double value) {
		// calculate sigma and generate 2 random variables
		double sigma = sqrt(1 / value);
		auto noise_out = noise.normal(sigma);

		// add noise 
		return { (in.first + noise_out.first) , (in.second + noise_out.second) };
	}

	// if hard decision then mapping back
	pair<int, int> hardDecision(pair<double, double>& in) {
		return { in.first >= 0 ? 0 : 1, in.second >= 0 ? 0 : 1 };
	}

	// encode one bit
	pair<pair<double, double>, int> encodeOneSoft(double SNR) {
		// get output of generator matrix
		auto out = generate(info_bits[ptr]);
		// exchange SNR value
		double value = pow(10, SNR / 10);
		// pass though AWGN channl
		auto out_AWGN = passAWGN(out, value);

		auto res = make_pair(out_AWGN, info_bits[ptr]);
		// recycle pointer
		ptr = (ptr + 1) % int(info_bits.size());
		return res;
	}

	pair<pair<int, int>, int> encodeOneHard(double SNR) {
		auto out = encodeOneSoft(SNR);
		return make_pair(hardDecision(out.first), out.second);
	}

	// encode all bits
	/*vector<pair<int, int>> encode(double SNR) {
		vector<pair<int, int>> res;

		for (int i = 0; i < N; i++) {
			if(i % 10000 == 0) cout << i << endl;
			res.push_back(encodeOne(SNR));
		}

		return res;
	}
	*/

	// output original info_bits to check decoded bits
	vector<int> outputOrigin(int num) {
		vector<int> res;
		for (int i = 0; i < num; i++) {
			int j = i % int(info_bits.size());
			res.push_back(info_bits[j]);
		}

		return res;
	}

	// iterate next new info_bit
	void genOneStep() {
		info_bits.push_back(*(info_bits.end() - 6) ^ *(info_bits.end() - 5));
	}

	// debugging function
	void checkGen(int num) {
		for (int i = 0; i < num; i++) {
			auto out = generate(info_bits[i]);
			cout << out.first << ' ' << out.second << endl;
		}
	}

	// debugging function
	void printInfo() {
		for (auto i : info_bits) cout << i << ' ';
		cout << endl;
	}
};

class Decoder {
	/**
		Since from any state e.g. 101010 we can only go to next 2 posible
		states by shift and add, say 110101 = (101010 >> 1) + (1 << 5)
		and 010101 = (101010 >> 1) + (0 << 5).
		So a 2D array(N * 2^6) is needed to backtrack the path
		graph array stored the minima error connection between nodes in
		2 layers(set -1 to mark unused states)

		For truncation(low space complexity) enable a double linked list
		for graph to reduce space complexity to O(32 * 2^6) = O(1)

		choose double linked list for backtrack path from end and low
		time complexity for pop_front and push_back operations O(1)

		dp array smallest cumulative number of error bits from
		first received signal to nth received signal

		for space complexity reduce dp from O(N^2) to O(N)

		param:
			nth : index of layers(for backtrack)
			states_num : 6
			total_states : 2^6
	*/

private:
	list<vector<int>> graph;
	vector<double> dp;
	int nth, states_num, total_states;

public:
	Decoder(int _states_num) : nth(0), states_num(_states_num) {
		total_states = 1 << states_num;
		graph.push_back(vector<int>(total_states, -1));
		dp = vector<double>(total_states);
	}

	bool isTop(int a, int b) {
		bitset<6> _a(a);
		bitset<6> _b(b);
		_reverse(_a);
		_reverse(_b);

		return _a.to_ulong() < _b.to_ulong();
	}

	void _reverse(bitset<6>& in) {
		for (int i = 0; i < 3; i++) {
			bool temp = in[i];
			in[i] = in[5 - i];
			in[5 - i] = temp;
		}
	}

	// given a received bit, build next layer of graph and calculate next dp array
	template <typename T>
	void forward(T& code, bool hard = false) {
		vector<double> _dp(total_states, INT_MAX);
		vector<int> next_layer(total_states, -1);

		auto Hamming = [&](pair<int, int>& in) {
			return abs(code.first - in.first) + abs(code.second - in.second);
		};

		auto mapping = [](int x) { return x == 0 ? 1 : -1; };

		auto Euclidean = [&](pair<int, int>& in) {
			in.first = mapping(in.first);
			in.second = mapping(in.second);

			return -(in.first * code.first + in.second * code.second);
		};

		// iterate all possible current states 
		for (int cur_state = 0; cur_state < total_states; cur_state++) {
			// check if it's used
			int last_state = graph.back()[cur_state];
			// for special case we need to go from 000000 state
			if (nth == 0) last_state = 0;
			// if not used there will not be some branch from unused states
			if (last_state == -1) continue;

			// change bits to array
			vector<int> states(6);
			for (int i = 0; i < 6; i++) {
				states[5 - i] = 1 & (cur_state >> i);
			}

			// there are 2 possible next states
			int next_1 = (cur_state >> 1) + (1 << (states_num - 1));
			int next_0 = cur_state >> 1;

			// generate corresponding outputs
			auto out_1 = generate(1, states);
			auto out_0 = generate(0, states);

			// calculate error bits number
			double error_1 = hard ? Hamming(out_1) : Euclidean(out_1);
			double error_0 = hard ? Hamming(out_0) : Euclidean(out_0);

			//_dp[next_0] = min(_dp[next_0], dp[cur_state] + error_0);
			//_dp[next_1] = min(_dp[next_1], dp[cur_state] + error_1);

			// transition function of DP and update the path with smallest error
			if (dp[cur_state] + error_0 < _dp[next_0])// || 
				//(dp[cur_state] + error_0 == _dp[next_0] && isTop(cur_state, next_layer[next_0]))) 
			{
				_dp[next_0] = dp[cur_state] + error_0;
				next_layer[next_0] = cur_state;
			}

			if (dp[cur_state] + error_1 < _dp[next_1]) //|| 
				//(dp[cur_state] + error_1 == _dp[next_1] && isTop(cur_state, next_layer[next_1]))) 
			{
				_dp[next_1] = dp[cur_state] + error_1;
				next_layer[next_1] = cur_state;
			}

			// special case ending
			if (nth == 0) break;
		}

		// build next layer
		graph.push_back(next_layer);
		nth++;
		// update dp array
		dp = _dp;
	}

	pair<int, int> generate(int bit, vector<int> states) {
		int x1 = bit ^ states[1] ^ states[2] ^ states[4] ^ states[5];
		int x2 = bit ^ states[0] ^ states[1] ^ states[2] ^ states[5];

		return { x1,  x2 };
	}

	// decode : backtrack from last layer to first layer and output path(decode bits)
	int decode() {
		// find best path end
		double minima_error = INT_MAX;
		int minima_state = -1;
		for (int i = 0; i < total_states; i++) {
			if (dp[i] < minima_error || (dp[i] == minima_error && isTop(i, minima_state))) {
				minima_error = dp[i];
				minima_state = i;
			}
		}

		int _nth = nth;
		auto it = graph.rbegin();

		// backtrack path layer by layer
		while (_nth > 1) {
			//res.push_back(1 & (minima_state >> 5));
			minima_state = (*it)[minima_state];
			it++;
			_nth--;
		}

		// need to reverse the backtracked path to be decode bits
		// reverse(res.begin(), res.end());

		// throw away decoded bit
		graph.pop_front();
		nth--;
		// we can inference info_bit which is the highest bit of next state
		return 1 & (minima_state >> (states_num - 1));
	}
};

int main() {
	while (1) {
		int N = 1000;
		int states_num = 6;
		int truncate = 31;
		double SNR_start = 0, SNR_end = 0;
		ULL SEED = 100000000000000;
		bool hard = false;
		char _in;

		cout << "Please give the number of decoded bits N : ";
		cin >> N;
		cout << "Please give SNR : from _ and to _ ";
		cin >> SNR_start >> SNR_end;
		cout << "Please give SEED : ";
		cin >> SEED;
		cout << "Hard?(Y/N) : ";
		cin >> _in;
		cout << "Truncate length: ";
		cin >> truncate;

		truncate--;
		hard = _in == 'Y' ? true : false;
		string filename = hard ? "Hard_res.txt" : "Soft_res.txt";
		ofstream file(filename);

		vector<double> SNRs;
		vector<double> BERs;
		vector<int> errors;

		for (double SNR = SNR_start; SNR <= SNR_end; SNR += 0.5) {
			SNRs.push_back(SNR);
			Encoder encoder(N + truncate, SEED);
			Decoder decoder(states_num);
			int cnt = 0;
			vector<int> res;
			list<int> encodes;

			for (int i = 0; i < N + truncate; i++) {
				if (hard) {
					auto code = encoder.encodeOneHard(SNR);
					encodes.push_back(code.second);

					decoder.forward<pair<int, int>>(code.first, true);
				}
				else {
					auto code = encoder.encodeOneSoft(SNR);
					encodes.push_back(code.second);

					decoder.forward<pair<double, double>>(code.first, false);
				}

				//if (i % 100000 == 0) cout << "Number of decoded bits : " << i << endl;

				if (i >= truncate) {
					res.push_back(decoder.decode());
					if (res.back() != *encodes.begin()) cnt++;
					encodes.pop_front();
				}
			}

			cout << "The number of error bits decoded : " << cnt << endl;
			cout << "BER : " << (double)cnt / N << endl;
			BERs.push_back((double)cnt / N);
			errors.push_back(cnt);
		}

		for (auto& SNR : SNRs) file << SNR << ' ';
		file << endl;
		for (auto& BER : BERs) file << BER << ' ';
		file << endl;
		for (auto& error : errors) file << error << ' ';
	}
	return 0;
}
