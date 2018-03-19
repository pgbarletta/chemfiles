%option reentrant
%option bison-bridge
%option noyywrap

%{

// $Id$
//
//  Copyright (C) 2001-2010 Randal Henne, Greg Landrum and Rational Discovery LLC
//
//  Copyright (C) 2018 Jonathan Fine (Chemfiles version)
//
//   @@ All Rights Reserved  @@
//

#include <cstdio>
#ifdef WIN32
#include <io.h>
#endif

#include <string>
#include <cstring>
#include <vector>
#include <list>
#include <map>

#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/ErrorFmt.hpp"

#include "smiles.tab.hpp"

using namespace chemfiles;

#define YY_FATAL_ERROR(msg) smiles_lexer_error(msg)

void smiles_lexer_error(const char *msg) {
     throw format_error(msg);
}

size_t setup_smiles_string(const std::string &text,yyscan_t yyscanner){
//  YY_BUFFER_STATE buff=yysmiles__scan_string(text.c_str()+pos,yyscanner);
  // Faster implementation of yysmiles__scan_string that handles trimming
  YY_BUFFER_STATE b;
  char *buf;
  yyconst char * yybytes = text.c_str();
  yy_size_t _yybytes_len=text.size(), n, start, end;
  /* Get memory for full buffer, including space for trailing EOB's. */
  n = _yybytes_len + 2;
  buf = reinterpret_cast<char *>(yysmiles_alloc(n ,yyscanner ));
  if ( ! buf )
    smiles_lexer_error( "out of dynamic memory in yysmiles__scan_bytes()" );

  // ltrim

  for(start = 0 ; start < _yybytes_len; ++start) {
    if (yybytes[start] > 32) break;
  }
  for(end = _yybytes_len ; end > start; --end) {
    if (yybytes[end] > 32) break;
  }

  _yybytes_len = end-start+1;
  n = _yybytes_len + 2;
  memcpy(buf, yybytes+start, _yybytes_len);


  buf[_yybytes_len] = buf[_yybytes_len+1] = YY_END_OF_BUFFER_CHAR;

  b = yysmiles__scan_buffer(buf,n ,yyscanner);
  if ( ! b )
    smiles_lexer_error( "bad buffer in yysmiles__scan_bytes()" );

  /* It's okay to grow etc. this buffer, and we should throw it
   * away when we're done.
   */
  b->yy_is_our_buffer = 1;


  if(!b)
    throw format_error("invalid buffer");
  return start;

}
%}

%s IN_ATOM_STATE
%%

@[' ']*TH |
@[' ']*AL |
@[' ']*SQ |
@[' ']*BP |
@[' ']*OH   { return CHI_CLASS_TOKEN; }

@           { return AT_TOKEN; }


<IN_ATOM_STATE>He       { yylval->atom = new Atom("", "He"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Li       { yylval->atom = new Atom("", "Li"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Be       { yylval->atom = new Atom("", "Be"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ne       { yylval->atom = new Atom("", "Ne"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Na       { yylval->atom = new Atom("", "Na"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mg       { yylval->atom = new Atom("", "Mg"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Al       { yylval->atom = new Atom("", "Al"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Si       { yylval->atom = new Atom("", "Si"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ar       { yylval->atom = new Atom("", "Ar"); return ATOM_TOKEN; }
<IN_ATOM_STATE>K        { yylval->atom = new Atom("", "K");  return ATOM_TOKEN; }
<IN_ATOM_STATE>Ca       { yylval->atom = new Atom("", "Ca"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sc       { yylval->atom = new Atom("", "Sc"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ti       { yylval->atom = new Atom("", "Ti"); return ATOM_TOKEN; }
<IN_ATOM_STATE>V        { yylval->atom = new Atom("", "V");  return ATOM_TOKEN; }
<IN_ATOM_STATE>Cr       { yylval->atom = new Atom("", "Cr"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mn       { yylval->atom = new Atom("", "Mn"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Fe       { yylval->atom = new Atom("", "Fe"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Co       { yylval->atom = new Atom("", "Co"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ni       { yylval->atom = new Atom("", "Ni"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cu       { yylval->atom = new Atom("", "Cu"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Zn       { yylval->atom = new Atom("", "Zn"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ga       { yylval->atom = new Atom("", "Ga"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ge       { yylval->atom = new Atom("", "Ge"); return ATOM_TOKEN; }
<IN_ATOM_STATE>As       { yylval->atom = new Atom("", "As"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Se       { yylval->atom = new Atom("", "Se"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Kr       { yylval->atom = new Atom("", "Kr"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rb       { yylval->atom = new Atom("", "Rb"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sr       { yylval->atom = new Atom("", "Sr"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Y        { yylval->atom = new Atom("", "Y");  return ATOM_TOKEN; }
<IN_ATOM_STATE>Zr       { yylval->atom = new Atom("", "Zr"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Nb       { yylval->atom = new Atom("", "Nb"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mo       { yylval->atom = new Atom("", "Mo"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Tc       { yylval->atom = new Atom("", "Tc"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ru       { yylval->atom = new Atom("", "Ru"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rh       { yylval->atom = new Atom("", "Rh"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pd       { yylval->atom = new Atom("", "Pd"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ag       { yylval->atom = new Atom("", "Ag"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cd       { yylval->atom = new Atom("", "Cd"); return ATOM_TOKEN; }
<IN_ATOM_STATE>In       { yylval->atom = new Atom("", "In"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sn       { yylval->atom = new Atom("", "Sn"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sb       { yylval->atom = new Atom("", "Sb"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Te       { yylval->atom = new Atom("", "Te"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Xe       { yylval->atom = new Atom("", "Xe"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cs       { yylval->atom = new Atom("", "Cs"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ba       { yylval->atom = new Atom("", "Ba"); return ATOM_TOKEN; }
<IN_ATOM_STATE>La       { yylval->atom = new Atom("", "La"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ce       { yylval->atom = new Atom("", "Ce"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pr       { yylval->atom = new Atom("", "Pr"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Nd       { yylval->atom = new Atom("", "Nd"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pm       { yylval->atom = new Atom("", "Pm"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sm       { yylval->atom = new Atom("", "Sm"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Eu       { yylval->atom = new Atom("", "Eu"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Gd       { yylval->atom = new Atom("", "Gd"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Tb       { yylval->atom = new Atom("", "Tb"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Dy       { yylval->atom = new Atom("", "Dy"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ho       { yylval->atom = new Atom("", "Ho"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Er       { yylval->atom = new Atom("", "Er"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Tm       { yylval->atom = new Atom("", "Tm"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Yb       { yylval->atom = new Atom("", "Yb"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Lu       { yylval->atom = new Atom("", "Lu"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Hf       { yylval->atom = new Atom("", "Hf"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ta       { yylval->atom = new Atom("", "Ta"); return ATOM_TOKEN; }
<IN_ATOM_STATE>W        { yylval->atom = new Atom("", "W");  return ATOM_TOKEN; }
<IN_ATOM_STATE>Re       { yylval->atom = new Atom("", "Re"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Os       { yylval->atom = new Atom("", "Os"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ir       { yylval->atom = new Atom("", "Ir"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pt       { yylval->atom = new Atom("", "Pt"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Au       { yylval->atom = new Atom("", "Au"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Hg       { yylval->atom = new Atom("", "Hg"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Tl       { yylval->atom = new Atom("", "Tl"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pb       { yylval->atom = new Atom("", "Pb"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Bi       { yylval->atom = new Atom("", "Bi"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Po       { yylval->atom = new Atom("", "Po"); return ATOM_TOKEN; }
<IN_ATOM_STATE>At       { yylval->atom = new Atom("", "At"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rn       { yylval->atom = new Atom("", "Rn"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Fr       { yylval->atom = new Atom("", "Fr"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ra       { yylval->atom = new Atom("", "Ra"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ac       { yylval->atom = new Atom("", "Ac"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Th       { yylval->atom = new Atom("", "Th"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pa       { yylval->atom = new Atom("", "Pa"); return ATOM_TOKEN; }
<IN_ATOM_STATE>U        { yylval->atom = new Atom("", "U");  return ATOM_TOKEN; }
<IN_ATOM_STATE>Np       { yylval->atom = new Atom("", "Np"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Pu       { yylval->atom = new Atom("", "Pu"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Am       { yylval->atom = new Atom("", "Am"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cm       { yylval->atom = new Atom("", "Cm"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Bk       { yylval->atom = new Atom("", "Bk"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cf       { yylval->atom = new Atom("", "Cf"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Es       { yylval->atom = new Atom("", "Es"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Fm       { yylval->atom = new Atom("", "Fm"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Md       { yylval->atom = new Atom("", "Md"); return ATOM_TOKEN; }
<IN_ATOM_STATE>No       { yylval->atom = new Atom("", "No"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Lr       { yylval->atom = new Atom("", "Lr"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rf       { yylval->atom = new Atom("", "Rf"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Db       { yylval->atom = new Atom("", "Db"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Sg       { yylval->atom = new Atom("", "Sg"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Bh       { yylval->atom = new Atom("", "Bh"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Hs       { yylval->atom = new Atom("", "Hs"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mt       { yylval->atom = new Atom("", "Mt"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ds       { yylval->atom = new Atom("", "Ds"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Rg       { yylval->atom = new Atom("", "Rg"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Cn       { yylval->atom = new Atom("", "Cn"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Nh       { yylval->atom = new Atom("", "Nh"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Fl       { yylval->atom = new Atom("", "Fl"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Mc       { yylval->atom = new Atom("", "Mc"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Lv       { yylval->atom = new Atom("", "Lv"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Ts       { yylval->atom = new Atom("", "Ts"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Og       { yylval->atom = new Atom("", "Og"); return ATOM_TOKEN; }

<IN_ATOM_STATE>Uun      { yylval->atom = new Atom("", "Uun"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uuu      { yylval->atom = new Atom("", "Uuu"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uub      { yylval->atom = new Atom("", "Uub"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uut      { yylval->atom = new Atom("", "Uut"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uuq      { yylval->atom = new Atom("", "Uuq"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uup      { yylval->atom = new Atom("", "Uup"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uuh      { yylval->atom = new Atom("", "Uuh"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uus      { yylval->atom = new Atom("", "Uus"); return ATOM_TOKEN; }
<IN_ATOM_STATE>Uuo      { yylval->atom = new Atom("", "Uuo"); return ATOM_TOKEN; }

B  { yylval->atom = new Atom("","B");return ORGANIC_ATOM_TOKEN; }
C  { yylval->atom = new Atom("","C");return ORGANIC_ATOM_TOKEN; }
N  { yylval->atom = new Atom("","N");return ORGANIC_ATOM_TOKEN; }
O  { yylval->atom = new Atom("","O");return ORGANIC_ATOM_TOKEN; }
P  { yylval->atom = new Atom("","P");return ORGANIC_ATOM_TOKEN; }
S  { yylval->atom = new Atom("","S");return ORGANIC_ATOM_TOKEN; }
F  { yylval->atom = new Atom("","F");return ORGANIC_ATOM_TOKEN; }
Cl { yylval->atom = new Atom("","Cl");return ORGANIC_ATOM_TOKEN; }
Br { yylval->atom = new Atom("","Br");return ORGANIC_ATOM_TOKEN; }
I  { yylval->atom = new Atom("","I");return ORGANIC_ATOM_TOKEN; }

H   {
        return H_TOKEN;
    }

b   {   yylval->atom = new Atom ( "", "B" );
        yylval->atom->set("aromatic", true);
        return AROMATIC_ATOM_TOKEN;
    }
c   {   yylval->atom = new Atom ( "", "C" );
        yylval->atom->set("aromatic", true);
        return AROMATIC_ATOM_TOKEN;
    }
n   {   yylval->atom = new Atom( "", "N" );
        yylval->atom->set("aromatic", true);
        return AROMATIC_ATOM_TOKEN;
    }
o   {   yylval->atom = new Atom( "", "O" );
        yylval->atom->set("aromatic", true);
        return AROMATIC_ATOM_TOKEN;
    }
p   {   yylval->atom = new Atom( "", "P" );
        yylval->atom->set("aromatic", true);
        return AROMATIC_ATOM_TOKEN;
    }
s   {   yylval->atom = new Atom( "", "S" );
        yylval->atom->set("aromatic", true);
        return AROMATIC_ATOM_TOKEN;
    }

<IN_ATOM_STATE>si   {   yylval->atom = new Atom( "", "Si" );
                        yylval->atom->set("aromatic", true);
                        return AROMATIC_ATOM_TOKEN;
                    }
<IN_ATOM_STATE>as   {   yylval->atom = new Atom( "", "As" );
                        yylval->atom->set("aromatic", true);
                        return AROMATIC_ATOM_TOKEN;
                    }
<IN_ATOM_STATE>se   {   yylval->atom = new Atom( "", "Se" );
                        yylval->atom->set("aromatic", true);
                        return AROMATIC_ATOM_TOKEN;
                    }
<IN_ATOM_STATE>te   {   yylval->atom = new Atom( "", "Te" );
                        yylval->atom->set("aromatic", true);
                        return AROMATIC_ATOM_TOKEN;
                    }

\*                  {   yylval->atom = new Atom();
                        yylval->atom->set("dummy", true);
                        // must be ORGANIC_ATOM_TOKEN because
                        // we aren't in square brackets:
                        return ORGANIC_ATOM_TOKEN;
                    }

<IN_ATOM_STATE>\:   { return COLON_TOKEN; }

%{
  // The next block is a workaround for a pathlogy in the SMILES produced
  // by some Biovia tools
%}
<IN_ATOM_STATE>\'Rf\'   { yylval->atom = new Atom("", "Rf"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Db\'   { yylval->atom = new Atom("", "Db"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Sg\'   { yylval->atom = new Atom("", "Sg"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Bh\'   { yylval->atom = new Atom("", "Bh"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Hs\'   { yylval->atom = new Atom("", "Hs"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Mt\'   { yylval->atom = new Atom("", "Mt"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Ds\'   { yylval->atom = new Atom("", "Ds"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Rg\'   { yylval->atom = new Atom("", "Rg"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Cn\'   { yylval->atom = new Atom("", "Cn"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Nh\'   { yylval->atom = new Atom("", "Nh"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Fl\'   { yylval->atom = new Atom("", "Fl"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Mc\'   { yylval->atom = new Atom("", "Mc"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Lv\'   { yylval->atom = new Atom("", "Lv"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Ts\'   { yylval->atom = new Atom("", "Ts"); return ATOM_TOKEN; }
<IN_ATOM_STATE>\'Og\'   { yylval->atom = new Atom("", "Og"); return ATOM_TOKEN; }

\=      {   yylval->bond = Bond::DOUBLE;
            return BOND_TOKEN; }
\#      {   yylval->bond = Bond::TRIPLE;
            return BOND_TOKEN; }
\:      {   yylval->bond = Bond::AROMATIC;
            //yylval->bond->setIsAromatic(true);
            return BOND_TOKEN; }
\-\>    {   yylval->bond = Bond::DATIVER;
            return BOND_TOKEN; }
\<\-    {   yylval->bond = Bond::DATIVEL;
            return BOND_TOKEN; }
\~      {   // This is a bond Query, not implemented
            return BOND_TOKEN;  }

[\\]{1,2}    { yylval->bond = Bond::DOWN;
        return BOND_TOKEN;  }

[\/]    { yylval->bond = Bond::UP;
        return BOND_TOKEN;  }

\-      { return MINUS_TOKEN; }

\+      { return PLUS_TOKEN; }

\(      { return GROUP_OPEN_TOKEN; }
\)      { return GROUP_CLOSE_TOKEN; }


\[      { BEGIN IN_ATOM_STATE; return ATOM_OPEN_TOKEN; }
<IN_ATOM_STATE>\]   { BEGIN INITIAL; return ATOM_CLOSE_TOKEN; }

\.      { return SEPARATOR_TOKEN; }

\%      { return PERCENT_TOKEN; }

[0]     { yylval->ival = 0; return ZERO_TOKEN; }
[1-9]   { yylval->ival = atoi( yytext ); return NONZERO_DIGIT_TOKEN; }



\n      return 0;

<<EOF>>     { return EOS_TOKEN; }
.             return yytext[0];

%%

#undef yysmiles_wrap
int yysmiles_wrap( void ) { return 1; }
