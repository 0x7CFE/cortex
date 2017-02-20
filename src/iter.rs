use std::iter::FromIterator;
use std::marker::PhantomData;

// fn main() {
//     for (i, chunk) in (1..8).collect_chunks::<Vec<_>>(3).enumerate() {
//         println!("{}: {:?}", i, chunk);
//     }
// }

pub trait CollectChunksExt: Iterator + Sized {
    fn collect_chunks<Chunk>(self, len: usize) -> CollectChunks<Chunk, Self>
        where Chunk: FromIterator<Self::Item>;
}

impl<Iter> CollectChunksExt for Iter where Iter: Iterator {
    fn collect_chunks<Chunk>(self, len: usize) -> CollectChunks<Chunk, Self>
        where Chunk: FromIterator<Self::Item>
    {
        CollectChunks {
            iter: self,
            len: len,
            marker: PhantomData,
        }
    }
}

pub struct CollectChunks<Chunk, Iter> {
    iter: Iter,
    len: usize,
    marker: PhantomData<Chunk>,
}

impl<Chunk, Iter> Iterator for CollectChunks<Chunk, Iter>
    where Chunk: FromIterator<Iter::Item>, Iter: Iterator
{
    type Item = Chunk;

    fn next(&mut self) -> Option<Self::Item> {
        let mut chunk = self.iter.by_ref().take(self.len).peekable();
        if chunk.peek().is_some() {
            Some(chunk.collect())
        } else {
            None
        }
    }
}
