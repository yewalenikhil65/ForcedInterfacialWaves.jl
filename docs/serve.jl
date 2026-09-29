import Pkg

const DOCS_DIR = @__DIR__

Pkg.activate(DOCS_DIR)
Pkg.instantiate()

using LiveServer
using Sockets

const PORT = 8001

let ip = try
        string(Sockets.getipaddr())
    catch
        nothing
    end
    println("="^60)
    println("Docs server starting on port $PORT")
    println("  Local:   http://localhost:$PORT")
    if ip !== nothing
        println("  Network: http://$ip:$PORT")
    else
        println("  Network: could not detect LAN IP automatically;")
        println("           run `ipconfig getifaddr en0` to find it.")
    end
    println("="^60)
end

servedocs(;
    foldername=DOCS_DIR,
    host="0.0.0.0",
    port=PORT,
)
