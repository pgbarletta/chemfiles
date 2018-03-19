%{

  // $Id$
  //
  //  Copyright (C) 2001-2016 Randal Henne, Greg Landrum and Rational Discovery LLC
  //
  //  Copyright (C) 2018 Jonathan Fine (Chemfiles version)
  //
  //   @@ All Rights Reserved  @@
  //

#include <cstring>
#include <vector>
#include <list>
#include <map>

#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/ErrorFmt.hpp"

#include "smiles.tab.hpp"

extern int yysmiles_lex(YYSTYPE *,void *);


using namespace chemfiles;
void
yysmiles_error( const char *input,
                chemfiles::Topology&,
                size_t&,
                std::vector<Residue>&,
                std::list<size_t>&,
                std::map<size_t, std::pair<size_t, chemfiles::Bond::Type>>&,
                void *,const char * msg )
{
  format_error("Error parsing: '{}' due to '{}'", input, msg);
}


%}

%define api.pure
%lex-param   {yyscan_t *scanner}
%parse-param {const char *input}
%parse-param {chemfiles::Topology& topol}
%parse-param {size_t& active_atom}
%parse-param {std::vector<chemfiles::Residue>& molList}
%parse-param {std::list<size_t>& branchPoints}
%parse-param {std::map<size_t, std::pair<size_t, chemfiles::Bond::Type>>& ringMarks}
%parse-param {void *scanner}

%union {
  int                      moli;
  chemfiles::Atom *        atom;
  chemfiles::Bond::Type    bond;
  int                      ival;
}

%token <atom> AROMATIC_ATOM_TOKEN ATOM_TOKEN ORGANIC_ATOM_TOKEN
%token <ival> NONZERO_DIGIT_TOKEN ZERO_TOKEN
%token GROUP_OPEN_TOKEN GROUP_CLOSE_TOKEN SEPARATOR_TOKEN LOOP_CONNECTOR_TOKEN
%token MINUS_TOKEN PLUS_TOKEN CHIRAL_MARKER_TOKEN CHI_CLASS_TOKEN CHI_CLASS_OH_TOKEN
%token H_TOKEN AT_TOKEN PERCENT_TOKEN COLON_TOKEN
%token <bond> BOND_TOKEN
%type <moli> cmpd mol
%type <atom> atomd element chiral_element h_element charge_element simple_atom
%type <ival>  nonzero_number number ring_number digit
%token ATOM_OPEN_TOKEN ATOM_CLOSE_TOKEN
%token EOS_TOKEN

%%

/* --------------------------------------------------------------- */
cmpd: mol
| cmpd error EOS_TOKEN{
  yyclearin;
  yyerrok;
  YYABORT;
}
| cmpd EOS_TOKEN {
  YYACCEPT;
}
| error EOS_TOKEN {
  yyclearin;
  yyerrok;
  YYABORT;
}
;

/* --------------------------------------------------------------- */
// FIX: mol MINUS DIGIT
mol: atomd {
  size_t sz = molList.size();
  molList.emplace_back( std::move(Residue("", molList.size())) );
  Residue& curMol = molList[ sz ];
  $1->set("SmilesStart",true);
  topol.add_atom(*$1);
  curMol.add_atom(topol.size() - 1);
  active_atom = topol.size() - 1;
  delete $1;
  $$ = sz;
}

| mol atomd       {
  Residue& mp = molList[$$];
  size_t atomIdx1 = active_atom;
  topol.add_atom(*$2);
  size_t atomIdx2 = topol.size() - 1;
  topol.add_bond(atomIdx1,atomIdx2);
  mp.add_atom(atomIdx2);
  active_atom = atomIdx2;
  delete $2;
}

| mol BOND_TOKEN atomd  {
  Residue& mp = molList[$$];
  size_t atomIdx1 = active_atom;
  topol.add_atom(*$3);
  size_t atomIdx2 = topol.size() - 1;
  topol.add_bond(atomIdx1,atomIdx2, $2);
  mp.add_atom(atomIdx2);
  active_atom = atomIdx2;
  delete $3;
}

| mol MINUS_TOKEN atomd {
  Residue& mp = molList[$$];
  size_t atomIdx1 = active_atom;
  topol.add_atom(*$3);
  size_t atomIdx2 = topol.size() - 1;
  topol.add_bond(atomIdx1,atomIdx2);
  mp.add_atom(atomIdx2);
  active_atom = atomIdx2;
  delete $3;
}

| mol SEPARATOR_TOKEN atomd {
  $3->set("SmilesStart",true);
  topol.add_atom(*$3);
  
  size_t sz = molList.size();
  molList.emplace_back(std::move(Residue("", molList.size())));
  Residue& curMol = molList[ sz ];
  curMol.add_atom(topol.size() - 1);
  active_atom = topol.size() - 1;
  delete $3;
}

| mol ring_number {
  
  if( ringMarks.find($2) == ringMarks.end() ) {
    ringMarks[$2] ={active_atom, Bond::UNDEFINED};
  } else {
    topol.add_bond(ringMarks[$2].first, active_atom, ringMarks[$2].second);
    ringMarks.erase($2);
  }
}

| mol BOND_TOKEN ring_number {

  if( ringMarks.find($3) == ringMarks.end() ) {
    ringMarks[$3] ={active_atom, $2};
  } else {
    topol.add_bond(ringMarks[$3].first, active_atom, ringMarks[$3].second);
    ringMarks.erase($3);
  }

}

| mol MINUS_TOKEN ring_number {

  if( ringMarks.find($3) == ringMarks.end() ) {
    ringMarks[$3] ={active_atom, Bond::SINGLE};
  } else {
    topol.add_bond(ringMarks[$3].first, active_atom, ringMarks[$3].second);
    ringMarks.erase($3);
  }
}

| mol GROUP_OPEN_TOKEN atomd {
  topol.add_atom(*$3);

  size_t atomIdx2 = topol.size() - 1;
  molList[static_cast<size_t>($$)].add_atom(atomIdx2);
  topol.add_bond(active_atom,atomIdx2);

  branchPoints.push_back(active_atom);
  active_atom = atomIdx2;
  delete $3;
}

| mol GROUP_OPEN_TOKEN BOND_TOKEN atomd  {
  topol.add_atom(*$4);

  size_t atomIdx2 = topol.size() - 1;
  molList[static_cast<size_t>($$)].add_atom(atomIdx2);
  topol.add_bond(active_atom,atomIdx2, $3);

  branchPoints.push_back(active_atom);
  active_atom = atomIdx2;
  delete $4;
}

| mol GROUP_OPEN_TOKEN MINUS_TOKEN atomd {
  topol.add_atom(*$4);

  size_t atomIdx2 = topol.size() - 1;
  molList[static_cast<size_t>($$)].add_atom(atomIdx2);
  topol.add_bond(active_atom,atomIdx2, Bond::SINGLE);

  branchPoints.push_back(active_atom);
  active_atom = atomIdx2;
  delete $4;
}

| mol GROUP_CLOSE_TOKEN {
  if(branchPoints.empty())
    yyerror(input,topol,active_atom,molList,branchPoints,ringMarks,scanner,"extra close parentheses");
  active_atom = branchPoints.back();
  branchPoints.pop_back();
}
;


/* --------------------------------------------------------------- */
atomd:	simple_atom
| ATOM_OPEN_TOKEN charge_element COLON_TOKEN number ATOM_CLOSE_TOKEN
{
  $$ = $2;
  $$->set("Implicit",false);
  $$->set("AtomMapNumber",$4);
}

| ATOM_OPEN_TOKEN charge_element ATOM_CLOSE_TOKEN
{
  $$ = $2;
  $2->set("Implicit",false);
}
;

/* --------------------------------------------------------------- */
charge_element:	h_element
| h_element PLUS_TOKEN { $1->set_charge(1); }
| h_element PLUS_TOKEN PLUS_TOKEN { $1->set_charge(2); }
| h_element PLUS_TOKEN number { $1->set_charge($3); }
| h_element MINUS_TOKEN { $1->set_charge(-1); }
| h_element MINUS_TOKEN MINUS_TOKEN { $1->set_charge(-2); }
| h_element MINUS_TOKEN number { $1->set_charge(-$3); }
		;

/* --------------------------------------------------------------- */
h_element:      H_TOKEN { $$ = new Atom("", "H"); }
                | number H_TOKEN { $$ = new Atom("", "H"); $$->set_mass($1); }
                | H_TOKEN H_TOKEN { $$ = new Atom("", "H"); $$->set("NumExplicitHs", 1); }
                | number H_TOKEN H_TOKEN { $$ = new Atom("", "H"); $$->set_mass($1); $$->set("NumExplicitHs", 1);}
                | H_TOKEN H_TOKEN number { $$ = new Atom("", "H"); $$->set("NumExplicitHs", $3); }
                | number H_TOKEN H_TOKEN number { $$ = new Atom("", "H"); $$->set_mass($1); $$->set("NumExplicitHs", $4);}
                | chiral_element
                | chiral_element H_TOKEN        { $$ = $1; $1->set("NumExplicitHs", 1);}
                | chiral_element H_TOKEN number	{ $$ = $1; $1->set("NumExplicitHs", $3);}
                ;

/* --------------------------------------------------------------- */
chiral_element:	 element
| element AT_TOKEN { $1->set("ChiralTag", "CHI_TETRAHEDRAL_CCW"); }
| element AT_TOKEN AT_TOKEN { $1->set("ChiralTag", "CHI_TETRAHEDRAL_CW"); }
;

/* --------------------------------------------------------------- */
element:	simple_atom
		|	number simple_atom { $2->set_mass( $1 ); $$ = $2; }
		|	ATOM_TOKEN
		|	number ATOM_TOKEN	   { $2->set_mass( $1 ); $$ = $2; }
		;

/* --------------------------------------------------------------- */
simple_atom:      ORGANIC_ATOM_TOKEN
                | AROMATIC_ATOM_TOKEN
                ;

/* --------------------------------------------------------------- */
ring_number:  digit
| PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit { $$ = $2*10+$3; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN { $$ = $3; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN { $$ = $3*10+$4; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*100+$4*10+$5; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*1000+$4*100+$5*10+$6; }
| PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*10000+$4*1000+$5*100+$6*10+$7; }
;

/* --------------------------------------------------------------- */
number:  ZERO_TOKEN
| nonzero_number
;

/* --------------------------------------------------------------- */
nonzero_number:  NONZERO_DIGIT_TOKEN
| nonzero_number digit { $$ = $1*10 + $2; }
;

digit: NONZERO_DIGIT_TOKEN
| ZERO_TOKEN
;

/*
  chival:	CHI_CLASS_TOKEN DIGIT_TOKEN
	| AT_TOKEN
        ;
*/


%%
