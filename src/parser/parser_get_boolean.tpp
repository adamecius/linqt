

template <typename T>
std::ifstream&
parser::GetBoolOption(std::ifstream& infile, const T text, bool& value)
{
	std::string Option;
	GetValue(infile,text,Option);
	std::transform(Option.begin(), Option.end(), Option.begin(), ::tolower);
	value = !((bool)Option.compare("1")*Option.compare("true"));
	return infile;
};
