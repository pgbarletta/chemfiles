// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cassert>
#include <sstream>
#include <list>
#include <map>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "chemfiles/formats/SMI.hpp"

#include "chemfiles/ErrorFmt.hpp"
#include "chemfiles/File.hpp"
#include "chemfiles/Frame.hpp"

int yysmiles_parse(const char *,
                   chemfiles::Topology&,
                   size_t&,
                   std::vector<chemfiles::Residue>&,
                   std::list<size_t>&,
                   std::map<size_t, std::pair<size_t, chemfiles::Bond::Type>>&,
                   void *);
int yysmiles_lex_init(void **);
int yysmiles_lex_destroy(void *);
size_t setup_smiles_string(const std::string &text, void *);
extern int yysmiles_debug;

using namespace chemfiles;

template<> FormatInfo chemfiles::format_information<SMIFormat>() {
    return FormatInfo("SMI").with_extension(".smi").description(
        "SMI text format"
    );
}

/// Fast-forward the file for one step, returning `false` if the file does
/// not contain one more step.
static bool forward(TextFile& file);

SMIFormat::SMIFormat(const std::string& path, File::Mode mode)
    : file_(TextFile::create(path, mode))
{
    while (!file_->eof()) {
        auto position = file_->tellg();
        if (!file_ || position == std::streampos(-1)) {
            throw format_error("IO error while reading '{}' as SMI", path);
        }
        if (forward(*file_)) {
            steps_positions_.push_back(position);
        }
    }
    file_->rewind();
}

size_t SMIFormat::nsteps() {
    return steps_positions_.size() - 1;
}

void SMIFormat::read_step(const size_t step, Frame& frame) {
    assert(step < steps_positions_.size());
    file_->seekg(steps_positions_[step]);
    read(frame);
}

void SMIFormat::read(Frame& frame) {

    Topology topol;
    size_t active_atom;
    std::vector<Residue> molVect;
    std::list<size_t> branchPoints;
    std::map<size_t, std::pair<size_t, Bond::Type>> ringMarks;
    void *scanner;

    if(yysmiles_lex_init(&scanner))
        throw memory_error("Error initializing lexer");
    try {
        auto line = file_->readline();
        size_t ltrim = setup_smiles_string(line, scanner);
        yysmiles_parse(line.c_str() + ltrim, topol, active_atom, molVect, branchPoints, ringMarks, scanner);
    } catch (...) {
        yysmiles_lex_destroy(scanner);
        throw;
    }
    yysmiles_lex_destroy(scanner);
    if (!branchPoints.empty()) {
        throw format_error("extra open parentheses");
    }

    for (auto res : molVect) {
        topol.add_residue(std::move(res));
    }

    frame.resize(topol.size());
    frame.set_topology(topol);
}

void SMIFormat::write(const Frame& frame) {
}

bool forward(TextFile& file) {
    if (!file) {return false;}

    try {
        auto line = file.readline();
    } catch (const FileError&) {
        // No more line left in the file
        return false;
    }
    
    return true;
}
