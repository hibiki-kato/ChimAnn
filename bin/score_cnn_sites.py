#!/usr/bin/env python3
"""Run Bowen's site scorer with metadata support for inference-only checkpoints.

Public checkpoints may contain epoch but omit optimizer_step. This adapter
records the unknown step as null, without changing weights or probabilities.
"""
from cnn_mamba import score_chromosome as scorer


def main():
    load, write = scorer.torch.load, scorer.atomic_json
    progress = {}

    def load_checkpoint(*args, **kwargs):
        checkpoint = load(*args, **kwargs)
        progress['epoch'] = checkpoint.get('epoch')
        checkpoint.setdefault('optimizer_step', None)
        return checkpoint

    def write_metadata(path, payload):
        if isinstance(payload, dict) and 'checkpoint_step' in payload:
            payload['checkpoint_epoch'] = progress.get('epoch')
        return write(path, payload)

    scorer.torch.load = load_checkpoint
    scorer.atomic_json = write_metadata
    try:
        scorer.main()
    finally:
        scorer.torch.load, scorer.atomic_json = load, write


if __name__ == '__main__':
    main()
