# 预编译负载只使用临时目录，避免污染包目录。
@setup_workload begin
    @compile_workload begin
        foldertest = joinpath(@__DIR__, "../test")
        # 预编译不得写入包目录：所有生成文件都放在临时目录，结束后清理
        tmpdir = mktempdir()
        try
            h, d = read_gadget2(joinpath(foldertest, "gassphere_littleendian.gadget2"), uAstro)
            write_gadget2(joinpath(tmpdir, "testGadget.gadget2"), h, d)
            write_gadget2(joinpath(tmpdir, "testGadgetHeaderGeneration.gadget2"), d)
            pos = read_gadget2_pos(joinpath(foldertest, "gassphere_littleendian.gadget2"), uAstro, uGadget2)
            h, d = read_gadget2(joinpath(foldertest, "pot_acc.format2.gadget2"), uAstro, uGadget2)
            h, d = read_gadget2(joinpath(foldertest, "pot_acc.format2.gadget2"), nothing, uGadget2)
            h, d = read_gadget2(joinpath(foldertest, "pot_acc.format2.gadget2"), nothing, nothing)
            pos = read_gadget2_pos(joinpath(foldertest, "pot_acc.format2.gadget2"), uGadget2)

            h, d = read_gadget2(joinpath(foldertest, "pot_acc.format2.gadget2"), uAstro, uGadget2)
            write_gadget2_format2(joinpath(tmpdir, "pot_acc.format2.test.gadget2"), h, d, acc = true, pot = true)
            # h, d = load(joinpath(foldertest, "gassphere_littleendian.gadget2"), uAstro)
            stars2d = [Star2D() for i = 1:10]
            write_csv(joinpath(tmpdir, "testcsvStar2D"), stars2d, nothing)
            stars = [Star(uAstro) for i = 1:10]
            write_csv(joinpath(tmpdir, "testcsvStar"), stars, uAstro)
            data = [[Star() for i = 1:10]; [Ball() for i = 1:10]]
            write_csv(joinpath(tmpdir, "testcsvGeneral"), data, nothing)
            write_ramses(joinpath(tmpdir, "ramses.csv"), stars)
            header, data = read_gadget2(joinpath(foldertest, "gassphere_littleendian.gadget2"), uAstro)
            d = [Star2D() for i = 1:10]
            write_gadget2_jld(joinpath(tmpdir, "testjldGadget.jld2"), header, d)
            h, d = read_gadget2_jld(joinpath(tmpdir, "testjldGadget.jld2"))
            write_jld(joinpath(tmpdir, "testjldGeneral.jld2"), d)
            d = read_jld(joinpath(tmpdir, "testjldGadget.jld2"))
        finally
            rm(tmpdir; recursive = true, force = true)
        end
    end
end
