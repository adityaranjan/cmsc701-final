#!/usr/bin/env bash
/root/.cargo/bin/cargo build --release
mv target/release/buildsa .
mv target/release/inspectsa .
mv target/release/querysa .
