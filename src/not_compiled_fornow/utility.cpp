




std::string printAttributes(std::unordered_map<std::string, std::string> a) {
	std::stringstream ss;
	for (auto it : a) {
		ss << it.first << "=" << it.second << ";";
	}
	return ss.str();
}

std::string printRecord(GFFRecord r) {
	std::string attr = printAttributes(r.attributes);

	std::stringstream ss;
	ss << "GFFRecord(seqid=" << r.seqid
		<< ", source=" << r.source 
		<< ", type=" << r.type
		<< ", start=" << r.start
		<< ", end=" << r.end 
		<< ", score=" << r.score 
		<< ", strand=" << r.strand
		<< ", phase=" << r.phase
		<< ", attributes={" << attr << "}";
	return ss.str();
}

template <typename T>
void printMap(std::unordered_map<std::string, T> &m, std::function<std::string (T)> f) {
	std::stringstream ss;
	for (auto it : m) {
		ss << it.first << ": " << f(it.second) << "\n";
	}
	Rcpp::Rcout << ss.str();
}

std::string VecToStr(std::vector<std::string> v) {
	std::stringstream ss;
	for (auto it : v) {
		ss << it << "\t";
	}
	return ss.str();
}

std::string PosToStr(Pos p) {
	std::stringstream ss;
	ss << "Pos(chr=" << p.chr;
	ss << ", start=" << p.start;
	ss << ", end=" << p.end;
	ss << ")";
	return ss.str();
}