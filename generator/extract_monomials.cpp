#include <mex.h>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <cctype>

void print_usage() {
	mexPrintf("A = extract_monomials(filename,n_eqs,n_basis,n_vars)\n");
}

class Monomial {
public:
	Monomial(int n_vars)
	{
		pow.resize(n_vars,0);
	}	
	void print() const {
		mexPrintf("[");
		for(int k = 0; k < pow.size(); k++) {
			mexPrintf("%d,",pow[k]);
		}
		mexPrintf("]\n");
	}
	std::vector<int> pow;
};



void debug_print(std::ifstream &input, std::string tag) {
	int pos0 = input.tellg();
	mexPrintf("DEBUG:(%s) '",tag.c_str());
	for(int k = 0; k < 20; k++) {		
		mexPrintf("%c",input.get());
	}
	mexPrintf("'\n");
	input.clear();
	input.seekg(pos0);	
}


int read_number(std::ifstream &input)
{
	int num;

	if(!isdigit(input.peek())) {
		mexPrintf("ERROR: Expected number. Actual: '%c'\n",input.peek());
		input.close();
		return 0;
	}


	input >> num;
	return num;
}

void parse_factor(std::ifstream &input, Monomial &output)
{
	// [\d* | x\d | x\d^\d ]
//	debug_print(input,"parse_factor");
	char c = input.peek();

	if(isdigit(c)) {
        if(c == '0') {
            // we have a zero element in the matrix
            input.get();
            output.pow.clear();
        } else {
            // read the rest of the numbers
            while(isdigit(input.peek()))
                input.get();	
        }
	} else if(c == 'x') {
		input.get();
		int k = read_number(input)-1;
		int p = 1;

		// check for higher powers
		if(input.peek() == '^') {
			input.get();
			p = read_number(input);
		}
		output.pow[k] = p;
	} else {
		mexPrintf("ERROR: expected number or 'x'. Actual: '%c\n'",c);
		input.close();
	}
}


void parse_term(std::ifstream &input, std::vector<Monomial> &output, int n_vars)
{
	// factor[\*factor]*
//	debug_print(input,"ENTER:parse_term");
	Monomial mon(n_vars);
	parse_factor(input,mon);
	while(input.peek() == '*') {
		input.get();
		parse_factor(input,mon);
	}
//	debug_print(input,"END:parse_term");
    if(mon.pow.size()>0)
        output.push_back(mon);
}
void parse_polynomial(std::ifstream &input, std::vector<Monomial> &output, int n_vars)
{
	// [-?]term([+|-]term)*
//	debug_print(input,"ENTER:parse_polynomial");
	// check for -
	if(input.peek() == '-')
		input.get();

	parse_term(input,output,n_vars);
	while(input.peek() == '+' || input.peek() == '-') {
		input.get();
		parse_term(input,output,n_vars);
	}
//	debug_print(input,"END:parse_polynomial");
}



void skip_whitespace(std::ifstream &input)
{
	while(isspace(input.peek()))
		input.get();
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs != 4 || nlhs != 1) {
		print_usage();
		return;
	}


	char* filename = mxArrayToString(prhs[0]);	
	int n_eqs = static_cast<int>(*mxGetPr(prhs[1]));
	int n_basis = static_cast<int>(*mxGetPr(prhs[2]));
	int n_vars = static_cast<int>(*mxGetPr(prhs[3]));


	char buf[8192];

	std::ifstream input(filename);
	input.rdbuf()->pubsetbuf(buf, 8192);

	input.seekg(9); // skip 'matrix {{'

	// create output array	
	mxArray* output_arr = mxCreateCellMatrix(n_eqs,n_basis);
	plhs[0] = output_arr;
	
	
	for(int eq_k = 0; eq_k < n_eqs; eq_k++) {		
		for(int basis_k = 0; basis_k < n_basis; basis_k++) {	

			std::vector<Monomial> mons;

			parse_polynomial(input,mons,n_vars);		

			// save result
			mxArray* m = mxCreateDoubleMatrix(n_vars, mons.size(), mxREAL);
			double *data = mxGetPr(m);
			int i = 0;
			for(int i = 0; i < mons.size(); i++) {
				for(int k = 0; k < n_vars; k++) {
					data[i*n_vars+k] = mons[i].pow[k];
				}				
  			}	  			
			int index = basis_k*n_eqs + eq_k;
			mxSetCell(output_arr,index,m);	


			if(input.peek() == ',') {
				input.get();
			}		
			skip_whitespace(input);
		}

		if(input.peek() != '}') {
			mexPrintf("ERROR: expected '}'. Actual '%c'\n",input.peek());
			input.close();
			return;
		}
		input.get();
		if(input.peek() == ',') {
			input.get(); // ','
			input.get(); // ' '	
			input.get(); // '{'	
		}		
	}

	
}
