// BSD 3-Clause License
//
// Copyright (c) 2023, Roy Ward
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// std::format support (C++20). A separate file, so <format> is only included when wanted

#pragma once

#ifndef PSEUDO_DOUBLE_FORMAT_CPP_H
#define PSEUDO_DOUBLE_FORMAT_CPP_H

#ifdef __has_include
#if __has_include(<version>)
#include <version>
#endif
#endif
#if !defined(__cpp_lib_format)
#error PseudoDouble_format.h needs C++20 std::format
#endif

#include "PseudoDouble.h"
#include <format>
#include <locale>
#include <string>
#include <string_view>

// The same format specs as double: [[fill]align][sign][#][0][width][.precision][L][type], with type e, E, f, F, g, G
// or none: none with no precision writes the shortest form that reads back (fixed or scientific, whichever is
// shorter, as std::format does for double). Not supported: a or A (hex), and width or precision given as nested
// arguments ({:{}}); both throw std::format_error. Up to 19 significant digits are correctly rounded (and integers
// below 2^127 written in fixed are exact); beyond that the digits are 0
template<> struct std::formatter<pseudodouble::PseudoDouble,char> {
	char fill=' ';
	char align=0;                                                      // '<', '>', '^', or 0 for the default (right)
	char sign='-';
	bool alt=false;
	bool zero=false;
	bool localized=false;
	int width=0;
	int precision=-1;
	char type=0;

	constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
		auto it=ctx.begin(), end=ctx.end();
		auto is_align=[](char c) {return c=='<' || c=='>' || c=='^';};
		if(it!=end && it+1!=end && is_align(it[1]) && *it!='{' && *it!='}') {
			fill=*it;
			align=it[1];
			it+=2;
		} else if(it!=end && is_align(*it)) {
			align=*it++;
		}
		if(it!=end && (*it=='+' || *it=='-' || *it==' ')) {
			sign=*it++;
		}
		if(it!=end && *it=='#') {
			alt=true;
			++it;
		}
		if(it!=end && *it=='0') {
			zero=true;
			++it;
		}
		while(it!=end && *it>='0' && *it<='9') {
			width=width*10+(*it++-'0');
		}
		if(it!=end && *it=='.') {
			++it;
			if(it==end || *it<'0' || *it>'9') {
				throw std::format_error("PseudoDouble: precision must be a number");
			}
			precision=0;
			while(it!=end && *it>='0' && *it<='9') {
				precision=precision*10+(*it++-'0');
			}
		}
		if(it!=end && *it=='L') {
			localized=true;
			++it;
		}
		if(it!=end && std::string_view("eEfFgG").find(*it)!=std::string_view::npos) {
			type=*it++;
		}
		if(it!=end && *it!='}') {
			throw std::format_error("PseudoDouble: unsupported format spec");
		}
		return it;
	}

	template<class FormatContext>
	typename FormatContext::iterator format(const pseudodouble::PseudoDouble& x, FormatContext& ctx) const {
		char conv=type==0?'g':type;
		int prec=precision;
		if(type!=0 && prec<0) {
			prec=6;                                                    // e, f and g default to 6, as for double
		}
		int flags=(sign=='+'?PD_FMT_PLUS:0)|(sign==' '?PD_FMT_SPACE:0)|(alt?PD_FMT_ALT:0);
		char point=localized?std::use_facet<std::numpunct<char>>(ctx.locale()).decimal_point():'.';
		char small[64];
		std::string s;
		int n=pdi_format_point(small,sizeof small,x.get_internal(),conv,prec,flags,point);
		if(n<(int)sizeof small) {
			s.assign(small,(size_t)n);
		} else {
			s.resize((size_t)n);
			pdi_format_point(&s[0],(size_t)n+1,x.get_internal(),conv,prec,flags,point);
		}
		if(width>(int)s.size()) {
			size_t pad=(size_t)width-s.size();
			bool nan=x.get_internal()==PD_NAN;
			if(zero && align==0 && !nan) {
				s.insert((s[0]=='-' || s[0]=='+' || s[0]==' ')?1:0,pad,'0'); // zeros go after the sign
			} else if(align=='<') {
				s.append(pad,fill);
			} else if(align=='^') {
				s.insert(0,pad/2,fill);
				s.append(pad-pad/2,fill);
			} else {
				s.insert(0,pad,fill);
			}
		}
		return std::copy(s.begin(),s.end(),ctx.out());
	}
};
#endif
