function Base.show(io::IO, h::HeaderGadget2)
    print(
        io, 
        "Gadget2 Header:\n",
        "  Particle Info :\n",
        "    |    Type   |  Amount, Mass\n",
        "    |-----------|-------------------------\n",
        "    |    Gas    | ", h.npart[1], ", ", h.mass[1] * 1.0e10u"Msun", "\n",
        "    |    Halo   | ", h.npart[2], ", ", h.mass[2] * 1.0e10u"Msun", "\n",
        "    |    Disk   | ", h.npart[3], ", ", h.mass[3] * 1.0e10u"Msun", "\n",
        "    |    Bulge  | ", h.npart[4], ", ", h.mass[4] * 1.0e10u"Msun", "\n",
        "    |    Star   | ", h.npart[5], ", ", h.mass[5] * 1.0e10u"Msun", "\n",
        "    | BlackHole | ", h.npart[6], ", ", h.mass[6] * 1.0e10u"Msun", "\n",
        "  Start time: ", h.time, "\n",
        "  Redshift: ", h.redshift, "\n",
    )
end