#ifndef _HON_H_
#define _HON_H_

#include <string>
#include <cctype>
#include <iostream>
#include <map>
#include <vector>
#include <cstdio>
#include <memory>
#include <fstream>

namespace hstd
{

namespace hon
{

	typedef std::string Key;

	enum ValueType
	{
		ValueTypeRealNumber,
		ValueTypeString,
		ValueTypeMap,
		ValueTypeArray,
		ValueTypeNull,
	};

	class IValue
	{
	private:
	public:
		virtual ValueType type() const = 0; 
	};
//	typedef IValue* Value;

	typedef std::shared_ptr<IValue> Value;

	class RealNumber : public IValue
	{
	private:
		const double value_;
	public:
		RealNumber() : value_(0.0) {}
		RealNumber(const double value) : value_(value) {}
		virtual ValueType type() const { return ValueType::ValueTypeRealNumber; }
		double value() const { return value_; }
	};

	class String : public IValue
	{
	private:
		const std::string str_;
	public:
		String() : str_("") {}
		String(const std::string& str) : str_(str) {}
		virtual ValueType type() const { return ValueType::ValueTypeString; }
		std::string string() const { return str_; }
	};

	class NullObject : public IValue
	{
	public:
		virtual ValueType type() const { return ValueType::ValueTypeNull; }
	};
	
	template <typename T>
	struct Response {
		bool valid;
		const T body;
		std::string info;

		Response() : valid(false), body(), info() {}
		Response(bool in_valid, T in_body, std::string &in_info) : valid(in_valid), body(in_body), info(in_info) {}

		Response<T> &operator=(const Response<T> &o) {
			valid = o.valid;
			body = o.body;
			info = o.info;
			return *this;
		}

		operator bool() const {
			return valid;
		}
		
		const T& operator*() {
			if (valid)
				return body;
			
			throw std::runtime_error("Invalid Response: " + info);
		}
	};
	
	template <>
	struct Response<RealNumber*> {
		bool valid;
		const RealNumber* body;
		std::string info;

		Response() : valid(false) {}
		Response(bool in_valid, RealNumber* in_body, std::string &in_info) : valid(in_valid), body(in_body), info(in_info) {}

		Response<RealNumber*> &operator=(const Response<RealNumber*> &o) {
			valid = o.valid;
			body = o.body;
			info = o.info;
			return *this;
		}

		operator bool() const {
			return valid;
		}

		const double operator*() {
			if (valid)
				return body->value();
			throw std::runtime_error("Invalid Realnumber: " + info);
		}
	};
	
	template <>
	struct Response<String*> {
		bool valid;
		const String* body;
		std::string info;

		Response() : valid(false) {}
		Response(bool in_valid, String* in_body, std::string &in_info) : valid(in_valid), body(in_body), info(in_info) {}

		Response<String*> &operator=(const Response<String*> &o) {
			valid = o.valid;
			body = o.body;
			info = o.info;
			return *this;
		}

		operator bool() const {
			return valid;
		}

		const std::string operator*() {
			if (valid)
				return body->string();
			throw std::runtime_error("Invalid String: " + info);
		}
	};
	
	typedef std::vector<Value> ValueArray;
	class Array : public IValue
	{
	private:
		ValueArray array_;
	public:
		virtual ValueType type() const { return ValueType::ValueTypeArray; }
		Array() {}
		Array(ValueArray& in_array) : array_(in_array) {}

		virtual ~Array() {}

		const ValueArray& valuearray() const {
			return array_;
		}

		Response<Value> at(const size_t index) {
			if (index >= array_.size()) {
				return Response<Value>(false, NULL, std::string(""));
			}

			return Response<Value>(true, array_[index], std::string(""));
		}
	};
	
	template <>
	struct Response<Array*> {
		bool valid;
		const Array* body;
		std::string info;

		Response() : valid(false) {}
		Response(bool in_valid, Array* in_body,  std::string &in_info) : valid(in_valid), body(in_body), info(in_info) {}

		Response<Array*> &operator=(const Response<Array*> &o) {
			valid = o.valid;
			body = o.body;
			info = o.info;
			return *this;
		}

		operator bool() const {
			return valid;
		}

		const ValueArray& operator*() const {
			if (valid)
				return body->valuearray();
			throw std::runtime_error("Invalid Array: " + info);
		}
	};

	typedef std::map<Key, Value> KeyValue;
	class Map : public IValue
	{
	private:
		KeyValue kv_map_;
	public:

		virtual ValueType type() const { return ValueType::ValueTypeMap; }
		
		Map() {}
		Map(KeyValue& kv_map) : kv_map_(kv_map) {}

		virtual ~Map() {}

		Response<Value> find(const Key& key) const {
			KeyValue::const_iterator ite = kv_map_.find(key);
			if (ite == kv_map_.end())
				return Response<Value>(false, NULL, std::string(key));

			return Response<Value>(true, ite->second, std::string(key));
		}

		Response<RealNumber*> realnumber(const Key& key) const {
			Response<Value> response = find(key);
			if (response.valid && response.body->type() == ValueTypeRealNumber) {
				return Response<RealNumber*>(true, dynamic_cast<RealNumber*>(response.body.get()), std::string(key));
			}
			return Response<RealNumber*>(false, NULL, std::string(key));
		}
		
		Response<String*> string(const Key& key) const {
			Response<Value> response = find(key);
			if (response.valid && response.body->type() == ValueTypeString) {
				return Response<String*>(true, dynamic_cast<String*>(response.body.get()), std::string(key));
			}
			return Response<String*>(false, NULL, std::string(key));
		}
		
		Response<Map*> map(const Key& key) const {
			Response<Value> response = find(key);
			if (response.valid && response.body->type() == ValueTypeMap) {
				return Response<Map*>(true, dynamic_cast<Map*>(response.body.get()), std::string(key));
			}
			return Response<Map*>(false, NULL, std::string(key));
		}
		
		Response<Array*> valuearray(const Key& key) const {
			Response<Value> response = find(key);
			if (response.valid && response.body->type() == ValueTypeArray) {
				return Response<Array*>(true, dynamic_cast<Array*>(response.body.get()), std::string(key));
			}
			return Response<Array*>(false, NULL, std::string(key));
		}

		std::string to_str(int indent = 0, int numspace = 4) const {
			std::string ret;

			std::string indent_str1;
			for (int i = 0; i < indent; ++i)
				indent_str1 += " ";
			std::string indent_str2;
			for (int i = 0; i < indent + numspace; ++i)
				indent_str2 += " ";

			ret += indent_str1 + "{\n";
			unsigned int counter = 0;
			for (KeyValue::const_iterator ite = kv_map_.begin(); ite != kv_map_.end(); ++ite) {

				ret += indent_str2 + "\"" + ite->first + "\": ";
				if (ite->second->type() == ValueTypeString) {
					ret += "\"" + dynamic_cast<String*>(ite->second.get())->string() + "\"";
				} else if (ite->second->type() == ValueTypeRealNumber) {
					ret += std::to_string((long double)dynamic_cast<RealNumber*>(ite->second.get())->value());
				} else if (ite->second->type() == ValueTypeMap) {
					ret += "\n" + dynamic_cast<Map*>(ite->second.get())->to_str(indent + numspace, numspace);
				}
				
				counter ++;
				if (counter < kv_map_.size())
					ret += ",";
				ret += "\n";
			}
			ret += indent_str1 + "}";

			return ret;
		}
	};

	
	typedef Response<Map*> ResponseMap;
	typedef Response<String*> ResponseString;
	typedef Response<RealNumber*> ResponseRealNumber;
	typedef Response<Array*> ResponseArray;
	
	typedef Response<Map> Root;
	class IParser
	{
	private:
	public:
		virtual Root parse(const std::string& input) = 0;
	};
	
	enum TokenType
	{
		TokenTypeString,
		TokenTypeRealNumber,
		TokenTypeComma,
		TokenTypeColon,
		TokenTypeBeginBrace,
		TokenTypeEndBrace,
		TokenTypeBeginBracket,
		TokenTypeEndBracket
	};

	class IToken
	{
	private:
	public:
		virtual enum TokenType type() const = 0;
		virtual std::string to_str() const = 0;
	};
	class TokenString : public IToken {
	private:
		const std::string str_;
	public: 
		TokenString(const std::string& str) : str_(str) {}
		virtual enum TokenType type() const { return TokenTypeString; } 
		std::string str() const { return str_; }
		virtual std::string to_str() const { return "TokenString(" + str_ + ")"; }
	};
	class TokenRealNumber : public IToken {
	private:
		const double value_;
	public: 
		TokenRealNumber(const double value) : value_(value) {}
		virtual enum TokenType type() const { return TokenTypeRealNumber; } 
		double value() const { return value_; }
		virtual std::string to_str() const { return "TokenRealNumber"; }
	};
	class TokenComma : public IToken {
	private:
	public: 
		virtual enum TokenType type() const { return TokenTypeComma; } 
		virtual std::string to_str() const { return "TokenComma"; }
	};
	class TokenColon : public IToken {
	private:
	public: 
		virtual enum TokenType type() const { return TokenTypeColon; } 
		virtual std::string to_str() const { return "TokenColon"; }
	};
	class TokenBeginBrace : public IToken {
	private:
	public: 
		virtual enum TokenType type() const { return TokenTypeBeginBrace; } 
		virtual std::string to_str() const { return "TokenBeginBrace"; }
	};
	class TokenEndBrace : public IToken {
	private:
	public: 
		virtual enum TokenType type() const { return TokenTypeEndBrace; } 
		virtual std::string to_str() const { return "TokenEndBrace"; }
	};
	class TokenBeginBracket : public IToken {
	private:
	public: 
		virtual enum TokenType type() const { return TokenTypeBeginBracket; } 
		virtual std::string to_str() const { return "TokenBeginBracket"; }
	};
	class TokenEndBracket : public IToken {
	private:
	public: 
		virtual enum TokenType type() const { return TokenTypeEndBracket; } 
		virtual std::string to_str() const { return "TokenEndBracket"; }
	};

	struct TokenInfo {
		IToken* token;
		// meta data
		unsigned int column;
		unsigned int row;

		TokenInfo(IToken *in_token, unsigned int in_column, unsigned int in_row) :
			token(in_token), column(in_column), row(in_row) {}
	};

	class StringParser : public IParser
	{
	private:
		bool error_output_;
		
		bool is_realnumber(const char ch) {
			return std::isdigit(ch) || ch == '.' ||  ch == '-' ||  ch == 'e';
		}

		Root parse_internal(const std::vector<TokenInfo>& tokens) {
			KeyValue tmp;
			const size_t tokens_length = tokens.size();
			TokenCursor cursor(tokens, 0, error_output_);

			tmp = block(&cursor);

			if (!cursor.valid() || !cursor.empty())
				return Root(false, Map(), std::string(""));

			return Root(true, tmp, std::string(""));
		}
		
		Response<std::vector<TokenInfo>> tokenize(const std::string& input) {
			std::vector<TokenInfo> ret;
			int now_index = 0;
			bool valid = true;
			const size_t input_length = input.size();

			int row = 0, column = 0;

			for (; now_index < input_length; ++now_index, ++row) {
				if (input[now_index] == '{') {
					ret.push_back(TokenInfo(new TokenBeginBrace, column, row));
				} else if (input[now_index] == '}') {
					ret.push_back(TokenInfo(new TokenEndBrace, column, row));
				} else if (input[now_index] == '[') {
					ret.push_back(TokenInfo(new TokenBeginBracket, column, row));
				} else if (input[now_index] == ']') {
					ret.push_back(TokenInfo(new TokenEndBracket, column, row));
				} else if (input[now_index] == ':') {
					ret.push_back(TokenInfo(new TokenColon, column, row));
				} else  if (input[now_index] == ',') {
					ret.push_back(TokenInfo(new TokenComma, column, row));
				} else if (input[now_index] == '"') {
					std::string tmp;
					now_index ++;
					row ++;
					for (; now_index < input_length; ++now_index, ++row) {
						if (input[now_index] == '"')
							break;
						tmp += input[now_index];
					}
					ret.push_back(TokenInfo(new TokenString(tmp), column, row));
				} else if (is_realnumber(input[now_index])) {
					std::string tmp;

					tmp += input[now_index];
					now_index ++;
					row ++;

					for (; now_index < input_length; ++now_index, ++row) {
						if (!is_realnumber(input[now_index])) {
							--now_index;
							--row;
							break;
						}
						tmp += input[now_index];
					}
					ret.push_back(TokenInfo(new TokenRealNumber(atof(tmp.c_str())), column, row));
				} else if (isspace(input[now_index])) {
					
					if (input[now_index] == '\n') {
						column ++;
						row = 0;
					}

				} else {
					valid = false;
					if (error_output_) {
						std::cerr << "Invalid token(" << column << "," << row << ") ";
						std::cerr << "\"" << input[now_index] << "\""  << std::endl;
					}
				}
			}

			return Response<std::vector<TokenInfo>>(valid, ret, std::string(""));
		}


		class TokenCursor
		{
		private:
			unsigned int now_index_;
			const std::vector<TokenInfo>& tokens_;
			bool error_;

			std::map<Key, Value> kv_;
			IToken *last_token_;

			bool error_output_;

		public:
			TokenCursor(const std::vector<TokenInfo>& tokens, const int now_index, bool error_output) :
				now_index_(now_index), tokens_(tokens), error_(false), last_token_(NULL), error_output_(error_output) {}

			IToken* last_token() {
				return last_token_;
			}

			bool empty() {
				return now_index_ >= tokens_.size();
			}

			bool valid() {
				return !error_;
			}

			bool accept(TokenType type) {
				if (now_index_ < tokens_.size() && tokens_[now_index_].token->type() == type) {
					// std::cout << tokens_[now_index_]->to_str() << std::endl;
					last_token_ = tokens_[now_index_].token;
					now_index_ ++;
					return true;
				}

				return false;
			}

			bool expect(TokenType type) {
				if (accept(type))
					return true;
				error_ = true;
				if (error_output_) {
					unsigned int row = 0;
					unsigned int column = 0;
					if (now_index_ < tokens_.size()) {
						row = tokens_[now_index_].row;
						column = tokens_[now_index_].column;
					} else if (tokens_.size() > 0) {
						row = tokens_[tokens_.size() - 1].row;
						column = tokens_[tokens_.size() - 1].column;
					}
					std::cerr << "Error(" << column << "," << row <<  "): unexpected token." << std::endl;
				}
				now_index_++;
				return false;
			}

			void append(Key& key, Value& value) {
			}
		};

		
		ValueArray bracket_body(TokenCursor *cursor) {
			ValueArray ret;

			for (;;) {
				// Value
				if (cursor->accept(TokenTypeRealNumber)) {
					const TokenRealNumber *token = dynamic_cast<TokenRealNumber*>(cursor->last_token());
					ret.push_back(Value(new RealNumber(token->value())));
				} else if (cursor->accept(TokenTypeString)) {
					const TokenString *token = dynamic_cast<TokenString*>(cursor->last_token());
					ret.push_back(Value(new String(token->str())));
				} else if (cursor->accept(TokenTypeBeginBracket)) {
					// [
					ValueArray tmp = bracket_body(cursor);

					// ]
					if (cursor->expect(TokenTypeEndBracket)) {
						ret.push_back(Value(new Array(tmp)));
					} else {
						break;
					}
				} else if (cursor->accept(TokenTypeBeginBrace)) {
					// {
					KeyValue tmp = brace_body(cursor);

					// }
					if (cursor->expect(TokenTypeEndBrace)) {
						ret.push_back(Value(new Map(tmp)));
					} else {
						break;
					}
				} else {
					break;
				}

				if (cursor->accept(TokenTypeComma)) {
				} else {
					break;
				}
			}

			return ret;
		}

		KeyValue brace_body(TokenCursor *cursor) {
			KeyValue ret;

			for (;;) {
				Key now_key;

				if (cursor->accept(TokenTypeString)) {
					// Key
					const TokenString *token = dynamic_cast<TokenString*>(cursor->last_token());
					now_key = token->str();
				} else {
					break;
				}

				if (cursor->expect(TokenTypeColon)) {
				} else {
					break;
				}
				
				// Value
				if (cursor->accept(TokenTypeRealNumber)) {
					const TokenRealNumber *token = dynamic_cast<TokenRealNumber*>(cursor->last_token());
					ret.insert(std::pair<Key, Value>(now_key, Value(new RealNumber(token->value()))));
				} else if (cursor->accept(TokenTypeString)) {
					const TokenString *token = dynamic_cast<TokenString*>(cursor->last_token());
					ret.insert(std::pair<Key, Value>(now_key, Value(new String(token->str()))));
				} else if (cursor->accept(TokenTypeBeginBracket)) {
					// [
					ValueArray tmp = bracket_body(cursor);

					// ]
					if (cursor->expect(TokenTypeEndBracket)) {
						ret.insert(std::pair<Key, Value>(now_key, Value(new Array(tmp))));
					} else {
						break;
					}
				} else if (cursor->expect(TokenTypeBeginBrace)) {
					// {
					KeyValue tmp = brace_body(cursor);

					// }
					if (cursor->expect(TokenTypeEndBrace)) {
						ret.insert(std::pair<Key, Value>(now_key, Value(new Map(tmp))));
					} else {
						break;
					}
				} else {
					break;
				} 

				if (cursor->accept(TokenTypeComma)) {
				} else {
					break;
				}
			}

			return ret;
		}

		KeyValue block(TokenCursor *cursor) {
			KeyValue ret;

			// {
			if (cursor->expect(TokenTypeBeginBrace)) {
			} else {
				return ret;
			}
			
			ret = brace_body(cursor);
			
			// }
			if (cursor->expect(TokenTypeEndBrace)) {
			} else {
				return ret;
			}
			return ret;
		}
	public:
		StringParser(bool error_output = true) : error_output_(error_output) {}
		virtual ~StringParser() {}

		virtual Root parse(const std::string& input)
		{
			std::string work;

			// Tokenize
			Response<std::vector<TokenInfo>> tokenize_response = tokenize(input);

			if (!tokenize_response.valid)
				return Root(false, Map(), std::string(""));

			std::vector<TokenInfo> tokens = tokenize_response.body;

			// Parse
			Root mp = parse_internal(tokens);
			
			for (unsigned int i = 0; i < tokens.size(); ++i) {
				delete tokens[i].token;
			}

			return mp;
		}
	};
	
	
	class FileParser : public IParser
	{
	private:
		bool error_output_;
		StringParser parser;
	public:
		FileParser(bool error_output = true) : error_output_(error_output), parser(error_output) {}
		
		virtual Root parse(const std::string& filename) {
			std::ifstream ifs(filename);
			if (ifs.fail()) {
				if (error_output_)
					std::cerr << "Error: file open failed(" << filename << ")" << std::endl;
				return Root(false, Map(), std::string(""));
			}

			std::string input;
			std::string line;
			while (std::getline(ifs, line)) {
				input += line + "\n";
			}

			return parser.parse(input);	
		}
	};

} // namespace hon

} // namespace hstd

#endif // _HON_H_