// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license
#include <catch.hpp>
#include <chemfiles.hpp>
using namespace chemfiles;

#undef assert
#define assert CHECK

TEST_CASE() {
    // [example]
    auto frame = Frame();
    frame.add_atom(Atom("Fe"), {0.0, 0.0, 0.0});
    frame.add_atom(Atom("Fe"), {1.0, 1.0, 1.0});
    frame.add_atom(Atom("Fe"), {2.0, 2.0, 2.0});

    for (Atom& atom: frame) {
        assert(atom.name() == "Fe");
    }
    // [example]
}
