FROM rust:1-bookworm AS builder

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release --locked --bin rustaroo

FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates cd-hit mafft \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/rustaroo /usr/local/bin/rustaroo
COPY LICENSE /usr/share/licenses/rustaroo/LICENSE

LABEL org.opencontainers.image.title="Rustaroo" \
      org.opencontainers.image.description="Rustaroo with CD-HIT and MAFFT" \
      org.opencontainers.image.source="https://github.com/Stgeorge2002/Rustaroo" \
      org.opencontainers.image.licenses="MIT"

WORKDIR /data
ENTRYPOINT ["rustaroo"]
CMD ["--help"]
