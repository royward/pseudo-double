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

#pragma once

#ifndef PSEUDO_DOUBLE_IOSTREAM_CPP_H
#define PSEUDO_DOUBLE_IOSTREAM_CPP_H

// Use a separate file for this, as don't want to always have iostream included to use this class

#include "PseudoDouble.h"
#include <ostream>
#include <locale>
#include <string>

// Writes x as a double would be written: fixed, scientific or default (%g) notation, precision, showpos, showpoint,
// uppercase, width, fill, left/right/internal, and the locale's decimal point. hexfloat is not supported, and writes
// the shortest form that reads back instead. Up to 19 significant digits are correctly rounded (and integers below
// 2^127 written in fixed are exact); beyond that the digits are 0
inline std::ostream &operator<<(std::ostream &os, pseudodouble::PseudoDouble const &x) {
	std::ios_base::fmtflags fl=os.flags();
	char conv='g';
	int precision=(int)os.precision();
	switch(fl&std::ios_base::floatfield) {
	case std::ios_base::fixed:
		conv='f';
		break;
	case std::ios_base::scientific:
		conv='e';
		break;
	case std::ios_base::fixed|std::ios_base::scientific:
		precision=-1;                                                  // hexfloat: not supported
		break;
	default:
		break;
	}
	if(fl&std::ios_base::uppercase) {
		conv=(char)(conv-'a'+'A');
	}
	int flags=((fl&std::ios_base::showpos)?PD_FMT_PLUS:0)|((fl&std::ios_base::showpoint)?PD_FMT_ALT:0);
	char point=std::use_facet<std::numpunct<char>>(os.getloc()).decimal_point();
	char small[64];
	std::string s;
	int n=pdi_format_point(small,sizeof small,x.get_internal(),conv,precision,flags,point);
	if(n<(int)sizeof small) {
		s.assign(small,(size_t)n);
	} else {
		s.resize((size_t)n);
		pdi_format_point(&s[0],(size_t)n+1,x.get_internal(),conv,precision,flags,point);
	}
	std::streamsize w=os.width();
	if(w>(std::streamsize)s.size()) {
		size_t pad=(size_t)w-s.size();
		char fill=os.fill();
		std::ios_base::fmtflags adj=fl&std::ios_base::adjustfield;
		if(adj==std::ios_base::left) {
			s.append(pad,fill);
		} else if(adj==std::ios_base::internal) {
			s.insert((s[0]=='-' || s[0]=='+')?1:0,pad,fill);           // padding goes after the sign
		} else {
			s.insert(0,pad,fill);
		}
	}
	os.width(0);
	return os<<s;
}
#endif
