use nvim_rs::{Neovim, compat::tokio::Compat, error::LoopError};
use tokio::io::{AsyncWrite, WriteHalf};
use tokio::task::JoinHandle;

pub struct NeovimConnection<T>
where
    T: AsyncWrite + Send + 'static,
{
    pub nvim: Neovim<Compat<WriteHalf<T>>>,
    pub io_handler: JoinHandle<Result<Result<(), Box<LoopError>>, tokio::task::JoinError>>,
    pub target: String,
}

impl<T> NeovimConnection<T>
where
    T: AsyncWrite + Send + 'static,
{
    pub fn new(
        nvim: Neovim<Compat<WriteHalf<T>>>,
        io_handler: JoinHandle<Result<Result<(), Box<LoopError>>, tokio::task::JoinError>>,
        target: String,
    ) -> Self {
        Self {
            nvim,
            io_handler,
            target,
        }
    }

    pub fn target(&self) -> &str {
        &self.target
    }
}
