#pragma once
#include "MACRO_DEF.h"
namespace pf{

	struct bool_elem {
		int index;
		bool value;
		bool_elem() {
			index = 0;
			value = 0;
		}
		bool_elem& operator=(const bool_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};

	class bool_box {
	public:
		std::vector<bool_elem> _bool_box;
		typedef std::vector<bool_elem>::iterator iterator;
		typedef std::vector<bool_elem>::const_iterator citerator;
		iterator  begin() { return _bool_box.begin(); };
		iterator  end() { return _bool_box.end(); };
		bool& operator[](const int index) {
			for (auto i = _bool_box.begin(); i < _bool_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			std::cout << "_bool_box error, can't find the value index : " << index << std::endl;
			SYS_PROGRAM_STOP;
		}
		bool_box& operator=(const bool_box& n) {
			_bool_box = n._bool_box;
			return *this;
		}
		void add_bool(int _index, bool _value) {
			for (auto i = _bool_box.begin(); i < _bool_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			_bool_box.reserve(_bool_box.size() + 1);
			bool_elem elem;
			elem.index = _index;
			elem.value = _value;
			_bool_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _bool_box.begin(); i < _bool_box.end();) {
				if (i->index == index) {
					i = _bool_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_bool_box.clear();
		}
		int size() const {
			return int(_bool_box.size());
		}

		bool_box() {
			_bool_box.reserve(0);
		}
		~bool_box() {
			_bool_box.clear();
		}
	};

	struct int_elem {
		int index;
		int value;
		int_elem() {
			index = 0;
			value = 0;
		}
		int_elem& operator=(const int_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};

	class int_box {
	public:
		int_box() {
			index = 0;
			_int_box.reserve(0);
		}
		~int_box() {
			_int_box.clear();
		}
		std::vector<int_elem> _int_box;
		int index;
		typedef std::vector<int_elem>::iterator iterator;
		typedef std::vector<int_elem>::const_iterator citerator;
		iterator  begin() { return _int_box.begin(); };
		iterator  end() { return _int_box.end(); };
		int& operator[](const int index) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			std::cout << "int_box error, can't find the value index : " << index << std::endl;
			SYS_PROGRAM_STOP;
		}
		int_box& operator=(const int_box& n) {
			_int_box = n._int_box;
			index = n.index;
			return *this;
		}
		void add_int(int _index, int _value) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			_int_box.reserve(_int_box.size() + 1);
			int_elem elem;
			elem.index = _index;
			elem.value = _value;
			_int_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _int_box.begin(); i < _int_box.end();) {
				if (i->index == index) {
					i = _int_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_int_box.clear();
		}
		int size() const {
			return int(_int_box.size());
		}
	};

	class int2_box {
	public:
		int2_box() {
			_int_box.reserve(0);
		}
		~int2_box() {
			_int_box.clear();
		}
		std::vector<int_box> _int_box;
		typedef std::vector<int_box>::iterator iterator;
		typedef std::vector<int_box>::const_iterator citerator;
		iterator  begin() { return _int_box.begin(); };
		iterator  end() { return _int_box.end(); };
		int_box& operator[](const int index) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i) {
				if (i->index == index) return (*i);
			}
			std::cout << "int2_box error, can't find the value index : " << index << std::endl;
			SYS_PROGRAM_STOP;
		}
		int2_box& operator=(const int2_box& n) {
			_int_box = n._int_box;
			return *this;
		}
		void add_int(int _index1, int _index2, int _value) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i)
				if (i->index == _index1) {
					i->add_int(_index2, _value);
					return;
				}
			_int_box.reserve(_int_box.size() + 1);
			int_box box;
			box.index = _index1;
			box.add_int(_index2, _value);
			_int_box.push_back(box);
		}
		void add_int(int _index) {
			for (auto i = _int_box.begin(); i < _int_box.end(); ++i)
				if (i->index == _index) {
					return;
				}
			_int_box.reserve(_int_box.size() + 1);
			int_box box;
			box.index = _index;
			_int_box.push_back(box);
		}
		void erase(int index) {
			for (auto i = _int_box.begin(); i < _int_box.end();) {
				if (i->index == index) {
					i = _int_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_int_box.clear();
		}
		int size() const {
			return int(_int_box.size());
		}
	};

	struct REAL_elem {
		int index;
		REAL value;
		REAL_elem() {
			index = 0;
			value = 0.0;
		}
		REAL_elem& operator=(const REAL_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};

	class REAL_box {
	public:
		REAL_box() {
			_REAL_box.reserve(0);
		}
		~REAL_box() {
			_REAL_box.clear();
		}
		std::vector<REAL_elem> _REAL_box;
		typedef std::vector<REAL_elem>::iterator iterator;
		typedef std::vector<REAL_elem>::const_iterator citerator;
		iterator  begin() { return _REAL_box.begin(); };
		iterator  end() { return _REAL_box.end(); };
		REAL& operator[](const int index) {
			for (auto i = _REAL_box.begin(); i < _REAL_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			std::cout << "REAL_box error, can't find the value index : " << index << std::endl;
			SYS_PROGRAM_STOP;
		}
		REAL_box& operator=(const REAL_box& n) {
			_REAL_box = n._REAL_box;
			return *this;
		}
		void add_REAL(int _index, REAL _value) {
			for (auto i = _REAL_box.begin(); i < _REAL_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			_REAL_box.reserve(_REAL_box.size() + 1);
			REAL_elem elem;
			elem.index = _index;
			elem.value = _value;
			_REAL_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _REAL_box.begin(); i < _REAL_box.end();) {
				if (i->index == index) {
					i = _REAL_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_REAL_box.clear();
		}
		int size() const {
			return int(_REAL_box.size());
		}
	};

	struct pair_flag {
		int index1;
		int index2;
		int value;
		pair_flag() {
			index1 = 0;
			index2 = 0;
			value = 0;
		}
		pair_flag& operator=(const pair_flag& n) {
			index1 = n.index1;
			index2 = n.index2;
			value = n.value;
			return *this;
		}
	};

	class pair_flag_box {
	public:
		pair_flag_box() {
			_flag_box.reserve(0);
		}
		~pair_flag_box() {
			_flag_box.clear();
		}
		std::vector<pair_flag> _flag_box;
		typedef std::vector<pair_flag>::iterator iterator;
		typedef std::vector<pair_flag>::const_iterator citerator;
		iterator  begin() { return _flag_box.begin(); };
		iterator  end() { return _flag_box.end(); };
		int& operator()(const int index1, const int index2) {
			for (auto i = _flag_box.begin(); i < _flag_box.end(); ++i) {
				if ((i->index1 == index1 && i->index2 == index2) || (i->index1 == index2 && i->index2 == index1)) return i->value;
			}
			std::cout << "_flag_box error, can't find the value indexs : " << index1 << ", " << index2 << std::endl;
			SYS_PROGRAM_STOP;
		}
		void operator=(const pair_flag_box& n) {
			_flag_box = n._flag_box;
		}
		void set_flag(int index1, int index2, int flag) {
			for (auto i = _flag_box.begin(); i < _flag_box.end(); ++i)
				if ((i->index1 == index1 && i->index2 == index2) || (i->index1 == index2 && i->index2 == index1)) {
					i->value = flag;
					return;
				}
			_flag_box.reserve(_flag_box.size() + 1);
			pair_flag elem;
			elem.index1 = index1;
			elem.index2 = index2;
			elem.value = flag;
			_flag_box.push_back(elem);
		}
		void erase(int index1, int index2) {
			for (auto i = _flag_box.begin(); i < _flag_box.end();) {
				if ((i->index1 == index1 && i->index2 == index2) || (i->index1 == index2 && i->index2 == index1)) {
					i = _flag_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_flag_box.clear();
		}
		int size() const {
			return int(_flag_box.size());
		}
	};

	struct string_elem {
		int index;
		std::string value;
		string_elem() {
			index = 0;
			value = "";
		}
		string_elem& operator=(const string_elem& n) {
			index = n.index;
			value = n.value;
			return *this;
		}
	};

	class string_box {
	public:
		string_box() {
			_string_box.reserve(0);
		}
		~string_box() {
			_string_box.clear();
		}
		std::vector<string_elem> _string_box;
		typedef std::vector<string_elem>::iterator iterator;
		typedef std::vector<string_elem>::const_iterator citerator;
		iterator  begin() { return _string_box.begin(); };
		iterator  end() { return _string_box.end(); };
		std::string& operator[](const int index) {
			for (auto i = _string_box.begin(); i < _string_box.end(); ++i) {
				if (i->index == index) return i->value;
			}
			std::cout << "string_box error, can't find the value index : " << index << std::endl;
			SYS_PROGRAM_STOP;
		}
		string_box& operator=(const string_box& n) {
			_string_box = n._string_box;
			return *this;
		}
		void add_string(int _index, std::string _value) {
			for (auto i = _string_box.begin(); i < _string_box.end(); ++i)
				if (i->index == _index) {
					i->value = _value;
					return;
				}
			_string_box.reserve(_string_box.size() + 1);
			string_elem elem;
			elem.index = _index;
			elem.value = _value;
			_string_box.push_back(elem);
		}
		void erase(int index) {
			for (auto i = _string_box.begin(); i < _string_box.end();) {
				if (i->index == index) {
					i = _string_box.erase(i);
				}
				else
					++i;
			}
		}
		void clear() {
			_string_box.clear();
		}
		int size() const {
			return int(_string_box.size());
		}
	};


}
