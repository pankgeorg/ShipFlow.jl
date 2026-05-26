#!/usr/bin/env julia
#
# I6: Post-process the F3 headline frames into a single animated GIF.
# Wraps the system `convert` (ImageMagick) since ffmpeg is not present
# in the workspace. Run this after `headline_demo.jl` has produced the
# per-frame PNGs.
#
# Env knobs:
#   FRAME_DIR  (default: runs/headline_demo/frames)
#   OUT_PATH   (default: runs/headline_demo/headline.gif)
#   DELAY      (default: 8 — centiseconds between frames; 8 ≈ 12.5 fps)
#   SCALE      (default: 50 — percent; smaller = smaller file)

using Printf

const FRAME_DIR = get(ENV, "FRAME_DIR",
    abspath(joinpath(@__DIR__, "..", "runs", "headline_demo", "frames")))
const OUT_PATH  = get(ENV, "OUT_PATH",
    abspath(joinpath(@__DIR__, "..", "runs", "headline_demo", "headline.gif")))
const DELAY = parse(Int, get(ENV, "DELAY", "8"))
const SCALE = parse(Int, get(ENV, "SCALE", "50"))

isdir(FRAME_DIR) || error("Frame dir not found: $FRAME_DIR")
frames = sort(filter(f -> endswith(f, ".png"),
              readdir(FRAME_DIR; join=true)))
isempty(frames) && error("No PNG frames in $FRAME_DIR")

# Sanity check ImageMagick.
try
    run(pipeline(`convert -version`, devnull))
catch
    error("`convert` not found on PATH. Install ImageMagick or write the " *
          "frames out using a different toolchain.")
end

@printf "Assembling %d frames → %s (delay=%d, scale=%d%%)\n" length(frames) OUT_PATH DELAY SCALE
cmd = `convert -delay $DELAY -loop 0 -resize $(SCALE)% $frames $OUT_PATH`
run(cmd)
bytes = stat(OUT_PATH).size
@printf "Done: %s  (%.2f MB)\n" OUT_PATH bytes/1e6
