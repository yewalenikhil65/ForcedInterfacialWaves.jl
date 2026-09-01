import Pkg

const DOCS_DIR = @__DIR__
const PACKAGE_ROOT = normpath(joinpath(DOCS_DIR, ".."))

# Use an absolute project path so Documenter's nested @example modules can
# resolve the local package even after Documenter changes directories.
Pkg.activate(DOCS_DIR)
Pkg.develop(Pkg.PackageSpec(path=PACKAGE_ROOT); io=devnull)
Pkg.instantiate()

using LiveServer

servedocs(;
    foldername=DOCS_DIR,
    include_dirs=[joinpath(PACKAGE_ROOT, "src")],
    host="0.0.0.0",
    port=8001,
)
